/*
 * This file is part of the TASCAR software, see <http://tascar.org/>
 *
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

#include "receivermod.h"
#include "hoa.h"

class hoa2d_hybrid_t : public TASCAR::receivermod_base_t {
public:
  class data_t : public TASCAR::receivermod_base_t::data_t {
  public:
    data_t(uint32_t order);
    // ambisonic weights:
    std::vector<float> B;
    int32_t prevK;
  };
  hoa2d_hybrid_t(tsccfg::node_t xmlsrc);
  ~hoa2d_hybrid_t();
  void add_pointsource(const TASCAR::pos_t& prel, double width, const TASCAR::wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t*);
  void add_diffuse_sound_field(const TASCAR::amb1wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t*);
  void configure();
  void postproc(std::vector<TASCAR::wave_t>& output);
  receivermod_base_t::data_t* create_state_data(double srate,uint32_t fragsize) const;
  int32_t order;
  uint32_t channels;
  HOA::encoder2D_t encode;
  std::vector<float> B;
  std::vector<float> deltaB;

  std::vector<int> vspkassign;
  float** vspkchannels; //array of pointers to float arrays
  std::vector<float> *lut; //array of float vectors
  int32_t diviser;
  float increment;
  double dt;
  uint32_t chunk_size;
};

hoa2d_hybrid_t::data_t::data_t(uint32_t channels )
{
  B = std::vector<float>(channels, 0.0f );
  prevK = 0;
}

hoa2d_hybrid_t::hoa2d_hybrid_t(tsccfg::node_t xmlsrc)
  : TASCAR::receivermod_base_t(xmlsrc),
  order(3)
{
  GET_ATTRIBUTE(order,"","Ambisonics order");

  if(order < 0)
    throw TASCAR::ErrMsg("Negative order is not possible.");
  encode.set_order(order);
  channels = 2*order + 1;
  B = std::vector<float>(channels, 0.0f);
  deltaB = std::vector<float>(channels, 0.0f);


  //Create lut of SH transforms
  GET_ATTRIBUTE(diviser,"","Diviser");
  if( diviser < 0 )
    throw TASCAR::ErrMsg("Negative diviser is not possible.");


  increment = 2*M_PI/diviser;
  //populate transform lut
  lut = new std::vector<float>[diviser];
  float az_inc = -M_PI;
  float el = 0.0f;
  for(int i = 0; i <diviser;++i){
    lut[i] = std::vector<float>(channels, 0.0f );
    encode(az_inc,el,lut[i]);
    az_inc += increment;
  }
}

hoa2d_hybrid_t::~hoa2d_hybrid_t()
{
  for (int i =0;i<diviser;i++){
    delete[] vspkchannels[i];
  }
  delete[] vspkchannels;
}

void hoa2d_hybrid_t::add_pointsource(const TASCAR::pos_t& prel, double width, const TASCAR::wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t* sd)
{
  data_t* state(dynamic_cast<data_t*>(sd));
  if(!state)
    throw TASCAR::ErrMsg("Invalid data type.");

  int32_t kmin(0);
  float az(prel.azim());
  kmin = (int)(0.5 + (az+M_PI)/increment);
  if (kmin == diviser){
    kmin = 0;
  }
  int32_t prevK = state->prevK;
  double sum = dt;
  if (kmin == prevK){
    vspkassign[kmin] =1;
    for( unsigned int i=0;i<chunk.size();i++){
      vspkchannels[kmin][i] += chunk[i];
    }
  }
  else{
    vspkassign[kmin] = 1;
    vspkassign[prevK] = 1;
    for( unsigned int i=0;i<chunk.size();i++){
      vspkchannels[kmin][i] += sum*chunk[i];
      vspkchannels[prevK][i] += (1-sum)*chunk[i];
      sum += dt;
    }
  }
  state->prevK = kmin;
}
void hoa2d_hybrid_t::configure()
{
  n_channels = channels;
  TASCAR::receivermod_base_t::configure();
  //dynamically create array to store output of each VSPK
  vspkchannels = new float*[diviser];
  for (int i = 0; i < diviser;i++){
    vspkchannels[i] = new float[n_fragment]{};
  }
  dt = 1.0/std::max( 1.0,(double )n_fragment);
  chunk_size = n_fragment;
  vspkassign = std::vector<int>(diviser, 0);
}
void hoa2d_hybrid_t::postproc(std::vector<TASCAR::wave_t>& output){
  uint32_t float_size = chunk_size*sizeof(*vspkchannels[0]);
  //encode each vspk channel into SH domain
  for (int i =0;i<diviser;i++){
    if (vspkassign[i] ==1){
      B = lut[i];
      for(uint32_t t=0;t<chunk_size;++t)
        for(uint32_t acn=0;acn<channels;++acn){
          output[acn][t] += B[acn]*vspkchannels[i][t];
        }
      memset(vspkchannels[i],0.0,float_size);
      vspkassign[i] = 0;
    }
  }

  TASCAR::receivermod_base_t::postproc(output);
}
void hoa2d_hybrid_t::add_diffuse_sound_field(const TASCAR::amb1wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t* sd)
{
  if( output.size() ){
    output[0].add( chunk.w(), sqrtf(2.0f) );
    if( output.size() > 3 ){
      output[1].add( chunk.y() );
      output[2].add( chunk.z() );
      output[3].add( chunk.x() );
    }
  }
}

TASCAR::receivermod_base_t::data_t* hoa2d_hybrid_t::create_state_data(double srate, uint32_t fragsize) const
{
  return new data_t(channels);
}

REGISTER_RECEIVERMOD(hoa2d_hybrid_t);

/*
 * Local Variables:
 * mode: c++
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * compile-command: "make -C .."
 * End:
 */
