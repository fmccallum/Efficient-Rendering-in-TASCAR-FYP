/*
 * This file is part of the TASCAR software, see <http://tascar.org/>
 *
 * Copyright (c) 2018 Giso Grimm
 * Copyright (c) 2019 Giso Grimm
 * Copyright (c) 2020 Giso Grimm
 * Copyright (c) 2021 Giso Grimm
 */
/*
 * TASCAR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, version 3 of the License.
 *
 * TASCAR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHATABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License, version 3 for more details.
 *
 * You should have received a copy of the GNU General Public License,
 * Version 3 along with TASCAR. If not, see <http://www.gnu.org/licenses/>.
 */

#include "errorhandling.h"
#include "scene.h"
#include "nmmintrin.h" // SSE4.2

#define LUT_SIZE 720

class nsplutSIMD_t : public TASCAR::receivermod_base_speaker_t {
public:
  class data_t : public TASCAR::receivermod_base_t::data_t {
  public:
    data_t( uint32_t chunksize,uint32_t channels );
    virtual ~data_t();

    uint32_t prevK;
    // point source speaker weights:
    float* point_w;
    float* point_dw;
    // ambisonic weights:
    float* diff_w;
    float* diff_dw;
    float* diff_x;
    float* diff_dx;
    float* diff_y;
    float* diff_dy;
    float* diff_z;
    float* diff_dz;
    double dt;
  };
  nsplutSIMD_t( tsccfg::node_t xmlsrc );
  virtual ~nsplutSIMD_t(){ }
  void add_pointsource( const TASCAR::pos_t& prel, double width, const TASCAR::wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t* );
  receivermod_base_t::data_t* create_state_data( double srate,uint32_t fragsize ) const;
  void add_variables( TASCAR::osc_server_t* srv );
  bool useall;

  
  int lut[LUT_SIZE]= {};
  double spacing = TASCAR_2PI/LUT_SIZE;


};


nsplutSIMD_t::data_t::data_t( uint32_t chunksize,uint32_t channels )
{
  point_w = new float[channels];
  point_dw = new float[channels];
  diff_w = new float[channels];
  diff_dw = new float[channels];
  diff_x = new float[channels];
  diff_dx = new float[channels];
  diff_y = new float[channels];
  diff_dy = new float[channels];
  diff_z = new float[channels];
  diff_dz = new float[channels];
  for( uint32_t k=0;k<channels;k++ )
    point_w[k] = point_dw[k] = diff_w[k] = diff_dw[k] = diff_x[k] = 
      diff_dx[k] = diff_y[k] = diff_dy[k] = diff_z[k] = diff_dz[k] = 0;
  prevK = 0;
  dt = 1.0/std::max( 1.0,(double )chunksize);
}

nsplutSIMD_t::data_t::~data_t()
{
  delete [] point_w;
  delete [] point_dw;
  delete [] diff_w;
  delete [] diff_dw;
  delete [] diff_x;
  delete [] diff_dx;
  delete [] diff_y;
  delete [] diff_dy;
  delete [] diff_z;
  delete [] diff_dz;
}

nsplutSIMD_t::nsplutSIMD_t( tsccfg::node_t xmlsrc )
  : TASCAR::receivermod_base_speaker_t( xmlsrc ),
  useall(false)
{
  GET_ATTRIBUTE_BOOL(useall,"activate all speakers independent of source position");
  // create a boost point list from speaker layout:
  
  double az = -TASCAR_PI;
  for (int i = 0;i<LUT_SIZE;i++){
    TASCAR::pos_t dir;
    dir.set_sphere(1,az,0);
    uint32_t kmin(0);
    double dmin( distance(dir,spkpos[kmin].unitvector ));
    double dist(0);
    for( unsigned int k=1;k<spkpos.size();k++)
      if( ( dist = distance(dir,spkpos[k].unitvector ))<dmin ){
        kmin = k;
        dmin = dist;
      }
    lut[i] = kmin;
    az+=spacing;
  }

  
}

void nsplutSIMD_t::add_variables( TASCAR::osc_server_t* srv )
{
  TASCAR::receivermod_base_speaker_t::add_variables( srv );
  srv->add_bool( "/useall", &useall );
}

void nsplutSIMD_t::add_pointsource( const TASCAR::pos_t& prel, double width, const TASCAR::wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t* sd )
{
  data_t* d( (data_t* )sd);
  TASCAR::pos_t psrc( prel.normal());

  float az(prel.azim());


  int index = static_cast<int>(0.5 + (TASCAR_PI+az)/spacing);
  if(index == LUT_SIZE){
    index = 0;
  }
  uint32_t kmin = lut[index];


  //find amount weights need to change on each sample of chunk to end at correct weights
  uint32_t prevK = d->prevK;
  if (kmin == prevK){
    unsigned int i = 0;
    for( ;i<chunk.size();i+=4){
      __m128 out_vec = _mm_loadu_ps(&output[kmin][i]);
      __m128 t_vec = _mm_loadu_ps(&chunk[i]);
      out_vec = _mm_add_ps(t_vec,out_vec);
      _mm_storeu_ps(&output[kmin][i], out_vec);

      //output[kmin][i] += chunk[i];
    }
    for(;i<chunk.size();i++){
      output[kmin][i] += chunk[i];
    }
  }
  else{
    unsigned int i = 0;
    __m128 dt_vec =_mm_set_ps(4*d->dt,4*d->dt,4*d->dt,4*d->dt);
    __m128 sum_vec =_mm_set_ps(4*d->dt,3*d->dt,2*d->dt,d->dt);
    __m128 one_vec =_mm_set_ps(1.0,1.0,1.0,1.0);
    for(;i<chunk.size();i+=4){
      __m128 out_vec = _mm_loadu_ps(&output[kmin][i]);
      __m128 pout_vec = _mm_loadu_ps(&output[prevK][i]);
      __m128 t_vec = _mm_loadu_ps(&chunk[i]);

      __m128 pt_vec = _mm_sub_ps(one_vec,sum_vec);
      pt_vec = _mm_mul_ps(pt_vec,t_vec);
      t_vec = _mm_mul_ps(t_vec,sum_vec);

      out_vec = _mm_add_ps(t_vec,out_vec);
      pout_vec = _mm_add_ps(pt_vec,pout_vec);
      
      _mm_storeu_ps(&output[kmin][i], out_vec);
      _mm_storeu_ps(&output[prevK][i], pout_vec);
      sum_vec = _mm_add_ps(sum_vec,dt_vec);
    }
    double sum = (i+1)*d->dt;
    for(;i<chunk.size();i++){
      output[kmin][i] += sum*chunk[i];
      output[prevK][i] += (1-sum)*chunk[i];
      sum += d->dt;
    }
  }
  d->prevK = kmin;


}

TASCAR::receivermod_base_t::data_t* nsplutSIMD_t::create_state_data( double srate,uint32_t fragsize ) const
{
  return new data_t( fragsize,spkpos.size() );
}

REGISTER_RECEIVERMOD( nsplutSIMD_t );

/*
 * Local Variables:
 * mode: c++
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * compile-command: "make -C .."
 * End:
 */
