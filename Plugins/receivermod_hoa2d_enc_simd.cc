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
#include "nmmintrin.h" // SSE4.2


class hoa2d_enc_simd_t : public TASCAR::receivermod_base_t {
public:
  class data_t : public TASCAR::receivermod_base_t::data_t {
  public:
    data_t(uint32_t order);
    // ambisonic weights:
    std::vector<float> B;
  };
  hoa2d_enc_simd_t(tsccfg::node_t xmlsrc);
  ~hoa2d_enc_simd_t();
  void add_pointsource(const TASCAR::pos_t& prel, double width, const TASCAR::wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t*);
  void add_diffuse_sound_field(const TASCAR::amb1wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t*);
  void configure() { n_channels = channels; };
  receivermod_base_t::data_t* create_state_data(double srate,uint32_t fragsize) const;
  int32_t order;
  uint32_t channels;
  HOA::encoder2D_t encode;
  std::vector<float> B;
  std::vector<float> deltaB;

};

hoa2d_enc_simd_t::data_t::data_t(uint32_t channels )
{
  B = std::vector<float>(channels, 0.0f );
}

hoa2d_enc_simd_t::hoa2d_enc_simd_t(tsccfg::node_t xmlsrc)
  : TASCAR::receivermod_base_t(xmlsrc),
  order(3)
{
  GET_ATTRIBUTE(order,"","Ambisonics order");
  if( order < 0 )
    throw TASCAR::ErrMsg("Negative order is not possible.");
  channels = 2*order+1;
  encode.set_order( order );
  B = std::vector<float>(channels, 0.0f );
  deltaB = std::vector<float>(channels, 0.0f );
}

hoa2d_enc_simd_t::~hoa2d_enc_simd_t()
{
}

void hoa2d_enc_simd_t::add_pointsource(const TASCAR::pos_t& prel, double width, const TASCAR::wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t* sd)
{

  data_t* state(dynamic_cast<data_t*>(sd));
  if( !state )
    throw TASCAR::ErrMsg("Invalid data type.");
  float az(prel.azim());
  float el(prel.elev());
  // calculate encoding weights:
  encode(az,el,B);
  // calculate incremental weights:

  const  __m128 inc_t = _mm_set1_ps(t_inc);
  uint32_t acn = 0;
  for(;acn+4<channels;acn+=4){
    __m128 b = _mm_loadu_ps(&B[acn]);
    __m128 sb = _mm_loadu_ps(&(state->B[acn]));
    b = _mm_sub_ps(b,sb);
    b = _mm_mul_ps(b,inc_t);
    _mm_storeu_ps(&deltaB[acn], b);
  }
  for(;acn<channels;acn++){
    deltaB[acn] = (B[acn] - state->B[acn])*t_inc;
  }

  for(uint32_t acn=0;acn<channels;++acn){
    float dB=deltaB[acn];
    //create b = [B[acn]+db + B + 2db B + 3db B + 4db]
    //create [5db 6db 7db 8db]
    __m128 b_vec = _mm_set_ps(B[acn]+4*dB,B[acn]+3*dB,B[acn]+2*dB,B[acn]+dB);
    __m128 db_vec =_mm_set_ps(7*dB,6*dB,5*dB,4*dB);
    uint32_t t=0;
    for(;t+4<chunk.size();t+=4){
      //get [t t+1 t+2 t+3]
      __m128 t_vec = _mm_loadu_ps(&chunk[t]);
      //calculate bxt
      t_vec = _mm_mul_ps(t_vec,b_vec);
      //get output[acn][t-t+4]
      __m128 out_vec = _mm_loadu_ps(&output[acn][t]);
      //output = output  + (bxt)
      out_vec = _mm_add_ps(t_vec,out_vec);
      _mm_storeu_ps(&output[acn][t], out_vec);
      //increment b 
      b_vec = _mm_add_ps(b_vec,db_vec);
    }
    state->B[acn] = state->B[acn]+ t*dB;
    for(;t<chunk.size();t++){
      output[acn][t] += (state->B[acn] += deltaB[acn]) * chunk[t];
    }
  }

  //for(uint32_t acn=0;acn<channels;++acn)
  //  state->B[acn] = B[acn];


}

void hoa2d_enc_simd_t::add_diffuse_sound_field(const TASCAR::amb1wave_t& chunk, std::vector<TASCAR::wave_t>& output, receivermod_base_t::data_t* sd)
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

TASCAR::receivermod_base_t::data_t* hoa2d_enc_simd_t::create_state_data(double srate, uint32_t fragsize) const
{
  return new data_t(channels);
}

REGISTER_RECEIVERMOD(hoa2d_enc_simd_t);

/*
 * Local Variables:
 * mode: c++
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * compile-command: "make -C .."
 * End:
 */
