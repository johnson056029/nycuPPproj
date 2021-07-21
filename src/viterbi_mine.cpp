//
//  viterbi.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "trellis.h"
#include <fstream>
#include <string>
// #include "nmmintrin.h"
#include <nmmintrin.h>
#include <chrono>

namespace StochHMM {
	
	void trellis::viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		//TODO: determine which model and chose the type of algorithm to use;
		viterbi();
	}
	
	void trellis::viterbi(){
		if (hmm->isBasic()){
			simple_viterbi();
		}
		else{
			fast_complex_viterbi();
		}
	}
	
	void trellis::simple_viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		simple_viterbi();
	}

        void trellis::trace_and_align( std::string seq_x, std::string seq_y, size_t tracer_x, size_t tracer_y, std::vector<std::vector<std::vector<float>>> &scoring_table_, std::vector<std::vector<int>> &traceback_table_, int seq_x_size_real){
//                 return;
//                 for(int i=0;i<=seq_x.size();i++){
//                         for(int j=0;j<=seq_x.size();j++){
//                                 std::cout << traceback_table_[i][j] << " ";
//                         }std::cout << std::endl;
//                 }std::cout << std::endl;
                // return;

                // size_t tracer_x = seq_x.size();
                // size_t tracer_x = seq_x_size_real;
                // size_t tracer_y = seq_y.size();
                enum enum_state {MMATCH, INST_X, INST_Y};
                std::vector<enum_state> path;
                // while( !( tracer_x == 1 && tracer_y == 1 )){
                while( tracer_x >= 1 && tracer_y >= 1 ){
                        if( traceback_table_[tracer_x][tracer_y] == MMATCH ){
                                tracer_x -= 1;
                                tracer_y -= 1;
                                path.push_back(MMATCH);
                        }
                        else if( traceback_table_[tracer_x][tracer_y] == INST_X ){
                                tracer_x -= 1;
                                path.push_back(INST_X);
                        }
                        else if( traceback_table_[tracer_x][tracer_y] == INST_Y ){
                                tracer_y -= 1;
                                path.push_back(INST_Y);
                        }
                }
                std::reverse( path.begin(), path.end() );

                std::vector<std::vector<char>> alignment (3, std::vector<char>());
                enum { SEQ_X, CONSS, SEQ_Y};
                tracer_x = tracer_y = 0;

                for(auto path_iter : path ){
                        if( path_iter == MMATCH ){
                                alignment[SEQ_X].push_back( seq_x[tracer_x] );
                                alignment[SEQ_Y].push_back( seq_y[tracer_y] );
                                alignment[CONSS].push_back( alignment[SEQ_X].back() == alignment[SEQ_Y].back() ? alignment[SEQ_X].back() : (scoring_table_[tracer_x][tracer_y][path_iter]>0 ? '+':' ') );
                                tracer_x += 1;
                                tracer_y += 1;
                        }
                        if( path_iter == INST_X ){
                                alignment[SEQ_X].push_back( seq_x[tracer_x] );
                                alignment[SEQ_Y].push_back('-');
                                alignment[CONSS].push_back(' '); 
                                tracer_x += 1;
                        }
                        if( path_iter == INST_Y ){
                                alignment[SEQ_X].push_back('-');
                                alignment[SEQ_Y].push_back( seq_y[tracer_y] );
                                alignment[CONSS].push_back(' '); 
                                tracer_y += 1;
                        }
                }


                for(auto iter1 : alignment){
                        for(auto iter2 : iter1)
                                std::cout << iter2;
                        std::cout << std::endl;
                }
                std::cout << std::endl;
        }

        void trellis::update_interval(int levels_tobe_updated, std::vector<std::vector<std::vector<float>>> &scoring_table_, std::vector<std::vector<int>> &traceback_table_, size_t iter_seq_y, int t, __m128 _d_4f, __m128 _e_4f){
                // std::cout << "i just want to write a function.\n";
                if( levels_tobe_updated < 1 )
                        return;

                int starting_level = 4 - levels_tobe_updated;
                int gap_index = starting_level * t + 1;

                float d_4f[4];
                float e_4f[4];
                _mm_store_ps( d_4f, _d_4f);
                _mm_store_ps( e_4f, _e_4f);

                // enum enum_state {MMATCH, INST_X, INST_Y};
                enum  {MMATCH, INST_X, INST_Y};

                // if( scoring_table_[gap_index][iter_seq_y][INST_X] >= scoring_table_[gap_index-1][iter_seq_y][MMATCH]-d_4f[0] &&
                if( scoring_table_[gap_index][iter_seq_y][INST_X] >= scoring_table_[gap_index-1][iter_seq_y][MMATCH]-d_4f[0] &&
                        scoring_table_[gap_index][iter_seq_y][INST_X] >= scoring_table_[gap_index-1][iter_seq_y][INST_X]-e_4f[0] ){

                        update_interval( -- levels_tobe_updated, scoring_table_, traceback_table_, iter_seq_y, t, _d_4f, _e_4f);
                }


                for(int i=1;i<=t;i++){

                        float score_INST_X_4f[4];
                        score_INST_X_4f[0] = scoring_table_[i+0*t][iter_seq_y][INST_X];
                        score_INST_X_4f[1] = scoring_table_[i+1*t][iter_seq_y][INST_X];
                        score_INST_X_4f[2] = scoring_table_[i+2*t][iter_seq_y][INST_X];
                        score_INST_X_4f[3] = scoring_table_[i+3*t][iter_seq_y][INST_X];

                        float up_score_MMATCH_4f[4];
                        up_score_MMATCH_4f[0] = scoring_table_[i+0*t-1][iter_seq_y][MMATCH];
                        up_score_MMATCH_4f[1] = scoring_table_[i+1*t-1][iter_seq_y][MMATCH];
                        up_score_MMATCH_4f[2] = scoring_table_[i+2*t-1][iter_seq_y][MMATCH];
                        up_score_MMATCH_4f[3] = scoring_table_[i+3*t-1][iter_seq_y][MMATCH];

                        float up_score_INST_X_4f[4];
                        up_score_INST_X_4f[0] = scoring_table_[i+0*t-1][iter_seq_y][INST_X];
                        up_score_INST_X_4f[1] = scoring_table_[i+1*t-1][iter_seq_y][INST_X];
                        up_score_INST_X_4f[2] = scoring_table_[i+2*t-1][iter_seq_y][INST_X];
                        up_score_INST_X_4f[3] = scoring_table_[i+3*t-1][iter_seq_y][INST_X];

                        __m128 _score_INST_X_4f = _mm_load_ps(score_INST_X_4f);
                        __m128 _up_score_MMATCH_4f = _mm_load_ps(up_score_MMATCH_4f);
                        __m128 _up_score_INST_X_4f = _mm_load_ps(up_score_INST_X_4f);

                        float tmp_4f[4];
                        __m128 _tmp_4f;

                // update from up case
                        __m128 _tmp1_4f = _mm_sub_ps(_up_score_MMATCH_4f, _d_4f);
                        __m128 _tmp2_4f = _mm_sub_ps(_up_score_INST_X_4f, _e_4f);
                        _tmp_4f = _mm_max_ps( _tmp1_4f, _tmp2_4f );
                        _mm_store_ps(tmp_4f, _tmp_4f); 

                        float max_4f[4];
                        __m128 _max_4f;

                        _max_4f = _mm_max_ps( _score_INST_X_4f, _tmp_4f );
                        _mm_store_ps( max_4f, _max_4f );

                        scoring_table_[i+0*t][iter_seq_y][INST_X] = max_4f[0];
                        scoring_table_[i+1*t][iter_seq_y][INST_X] = max_4f[1];
                        scoring_table_[i+2*t][iter_seq_y][INST_X] = max_4f[2];
                        scoring_table_[i+3*t][iter_seq_y][INST_X] = max_4f[3];

//                         scoring_table_[i+0*t][iter_seq_y][INST_X] = std::max( {scoring_table_[i+0*t][iter_seq_y][INST_X], tmp_4f[0] });
//                         scoring_table_[i+1*t][iter_seq_y][INST_X] = std::max( {scoring_table_[i+1*t][iter_seq_y][INST_X], tmp_4f[1] });
//                         scoring_table_[i+2*t][iter_seq_y][INST_X] = std::max( {scoring_table_[i+2*t][iter_seq_y][INST_X], tmp_4f[2] });
//                         scoring_table_[i+3*t][iter_seq_y][INST_X] = std::max( {scoring_table_[i+3*t][iter_seq_y][INST_X], tmp_4f[3] });

                }

                update_interval( -- levels_tobe_updated, scoring_table_, traceback_table_, iter_seq_y, t, _d_4f, _e_4f);
                
                return;
        }
	
	void trellis::simple_viterbi(){

                std::ifstream ifs_x ( "./examples/HBA_HUMAN.fa", std::ifstream::in);
                std::string tmpstr;
                std::string seq_x;
                while(ifs_x>>tmpstr)
                        seq_x.append(tmpstr);
                int seq_x_size_real = seq_x.size();
                int t = seq_x.size() / 4;
                t = seq_x.size() % 4 == 0 ? t : t+1;
                for(int i=0; i<seq_x.size()%4; i++)
                        seq_x.append("A");

                std::ifstream ifs_y ( "./examples/LGB2_LUPLU.fa", std::ifstream::in);
                std::string seq_y;
                while(ifs_y>>tmpstr)
                        seq_y.append(tmpstr);

                std::vector<std::vector<float>> BLOSUM62 = hmm->log_odds_matrix("BLOSUM62");

                enum  AminoAcid {A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V};

                std::map<char, int> amap
                {   
                        {'A',A},{'R',R},{'N',N},{'D',D},{'C',C},{'Q',Q},{'E',E},{'G',G},{'H',H},{'I',I},{'L',L},{'K',K},{'M',M},{'F',F},{'P',P},{'S',S},{'T',T},{'W',W},{'Y',Y},{'V',V},
                };

                float delta =   0.2;
                float epsilon = 0.4;
                float tau =     0.1;
                float eta =     0.5;

                std::vector<std::vector<float>> substitution_matrix (20, std::vector<float>(20, 0));
                for(int row_iter=0; row_iter<20; row_iter++ ){
                        for(int col_iter=0; col_iter<20; col_iter++ ){
                                substitution_matrix[row_iter][col_iter] = BLOSUM62[row_iter][col_iter] + log( (1-2*delta-tau)/(pow(1-eta,2)) );
                        }
                }
                
                float d = -log( (delta*(1-epsilon-tau))/((1-eta)*(1-2*delta-tau)) );
                float e = -log( (epsilon)/(1-tau) );
                float d_4f[4];
                float e_4f[4];
                for(int i=0;i<4;i++){
                        d_4f[i] = d;
                        e_4f[i] = e;
                }
                __m128 _d_4f = _mm_load_ps(d_4f);
                __m128 _e_4f = _mm_load_ps(e_4f);


                // std::vector<std::vector<std::vector<float>>> scoring_table_ (seq_x.size()+1, std::vector<std::vector<float>> (seq_y.size()+1, std::vector<float>(3, -INFINITY)));
                std::vector<std::vector<std::vector<float>>> scoring_table_ (4*t+1, std::vector<std::vector<float>> (seq_y.size()+1, std::vector<float>(3, -INFINITY)));
                enum enum_state {MMATCH, INST_X, INST_Y};
                std::vector<std::vector<int>> traceback_table_ (4*t+1, std::vector<int> (seq_y.size()+1, -1));

                scoring_table_[0][0][MMATCH] = 2 * log(eta);
                // for(size_t iter_seq_x=1; iter_seq_x < seq_x.size()+1; iter_seq_x++){


// TIMING START
                auto __begin__ = std::chrono::high_resolution_clock::now();


                for(size_t iter_seq_y=1; iter_seq_y <= seq_y.size(); iter_seq_y++){
                        for(size_t iter_seq_x=1; iter_seq_x <= t; iter_seq_x++){

                                // substitution score
                                float subscore_4f[4];
                                subscore_4f[0] = substitution_matrix[amap[seq_x[iter_seq_x+0*t]]][amap[seq_y[iter_seq_y]]];
                                subscore_4f[1] = substitution_matrix[amap[seq_x[iter_seq_x+1*t]]][amap[seq_y[iter_seq_y]]];
                                subscore_4f[2] = substitution_matrix[amap[seq_x[iter_seq_x+2*t]]][amap[seq_y[iter_seq_y]]];
                                subscore_4f[3] = substitution_matrix[amap[seq_x[iter_seq_x+3*t]]][amap[seq_y[iter_seq_y]]];

                                __m128 _subscore_4f = _mm_load_ps(subscore_4f);

                                // scores blocks
                                float score_4f[3][4];
                                float left_score_4f[3][4];
                                float up_score_4f[3][4];
                                float upleft_score_4f[3][4];

                                __m128 _score_4f[3];
                                __m128 _left_score_4f[3];
                                __m128 _up_score_4f[3];
                                __m128 _upleft_score_4f[3];

                                for(int iter_state=0; iter_state<3; iter_state++){
                                        // score
                                        score_4f[iter_state][0] = scoring_table_[iter_seq_x+0*t][iter_seq_y][iter_state];
                                        score_4f[iter_state][1] = scoring_table_[iter_seq_x+1*t][iter_seq_y][iter_state];
                                        score_4f[iter_state][2] = scoring_table_[iter_seq_x+2*t][iter_seq_y][iter_state];
                                        score_4f[iter_state][3] = scoring_table_[iter_seq_x+3*t][iter_seq_y][iter_state];
                                        // left score
                                        left_score_4f[iter_state][0] = scoring_table_[iter_seq_x+0*t][iter_seq_y-1][iter_state];
                                        left_score_4f[iter_state][1] = scoring_table_[iter_seq_x+1*t][iter_seq_y-1][iter_state];
                                        left_score_4f[iter_state][2] = scoring_table_[iter_seq_x+2*t][iter_seq_y-1][iter_state];
                                        left_score_4f[iter_state][3] = scoring_table_[iter_seq_x+3*t][iter_seq_y-1][iter_state];
                                        // upper score
                                        up_score_4f[iter_state][0] = scoring_table_[iter_seq_x+0*t-1][iter_seq_y][iter_state];
                                        up_score_4f[iter_state][1] = scoring_table_[iter_seq_x+1*t-1][iter_seq_y][iter_state];
                                        up_score_4f[iter_state][2] = scoring_table_[iter_seq_x+2*t-1][iter_seq_y][iter_state];
                                        up_score_4f[iter_state][3] = scoring_table_[iter_seq_x+3*t-1][iter_seq_y][iter_state];
                                        // upper left score
                                        upleft_score_4f[iter_state][0] = scoring_table_[iter_seq_x+0*t-1][iter_seq_y-1][iter_state];
                                        upleft_score_4f[iter_state][1] = scoring_table_[iter_seq_x+1*t-1][iter_seq_y-1][iter_state];
                                        upleft_score_4f[iter_state][2] = scoring_table_[iter_seq_x+2*t-1][iter_seq_y-1][iter_state];
                                        upleft_score_4f[iter_state][3] = scoring_table_[iter_seq_x+3*t-1][iter_seq_y-1][iter_state];

                                        // load the above values into the simd refisters
                                        _score_4f[iter_state] = _mm_load_ps(score_4f[iter_state]);
                                        _left_score_4f[iter_state] = _mm_load_ps(left_score_4f[iter_state]);
                                        _up_score_4f[iter_state] = _mm_load_ps(up_score_4f[iter_state]);
                                        _upleft_score_4f[iter_state] = _mm_load_ps(upleft_score_4f[iter_state]);
                                }

                                float tmp_4f[4];
                                __m128 _tmp_4f;

                                float max_4f[4];
                                __m128 _max_4f;


                                // for MATCH/MISMATCH case
                                // update by subsitution case
                                // prepare the maximum value plus substitution score for the upperleft entry
                                float max_upleft_score_4f[4];
                                __m128 _max_upleft_score_4f = _mm_max_ps( _upleft_score_4f[MMATCH], _upleft_score_4f[INST_X]);
                                _max_upleft_score_4f = _mm_max_ps( _max_upleft_score_4f, _upleft_score_4f[INST_Y]);
                                _tmp_4f = _mm_add_ps( _max_upleft_score_4f, _subscore_4f);
                                _mm_store_ps(tmp_4f, _tmp_4f); 

                                _max_4f = _mm_max_ps(_score_4f[MMATCH], _tmp_4f);
                                _mm_store_ps(max_4f, _max_4f);

                                scoring_table_[iter_seq_x+0*t][iter_seq_y][MMATCH] = max_4f[0]; 
                                scoring_table_[iter_seq_x+1*t][iter_seq_y][MMATCH] = max_4f[1]; 
                                scoring_table_[iter_seq_x+2*t][iter_seq_y][MMATCH] = max_4f[2]; 
                                scoring_table_[iter_seq_x+3*t][iter_seq_y][MMATCH] = max_4f[3]; 

//                                 scoring_table_[iter_seq_x+0*t][iter_seq_y][MMATCH] = std::max({ scoring_table_[iter_seq_x+0*t][iter_seq_y][MMATCH], tmp_4f[0] });
//                                 scoring_table_[iter_seq_x+1*t][iter_seq_y][MMATCH] = std::max({ scoring_table_[iter_seq_x+1*t][iter_seq_y][MMATCH], tmp_4f[1] });
//                                 scoring_table_[iter_seq_x+2*t][iter_seq_y][MMATCH] = std::max({ scoring_table_[iter_seq_x+2*t][iter_seq_y][MMATCH], tmp_4f[2] });
//                                 scoring_table_[iter_seq_x+3*t][iter_seq_y][MMATCH] = std::max({ scoring_table_[iter_seq_x+3*t][iter_seq_y][MMATCH], tmp_4f[3] });


                                // update from left case
                                __m128 _tmp1_4f;
                                __m128 _tmp2_4f;

                                _tmp1_4f = _mm_sub_ps(_left_score_4f[MMATCH], _d_4f);
                                _tmp2_4f = _mm_sub_ps(_left_score_4f[INST_Y], _e_4f);
                                _tmp_4f = _mm_max_ps( _tmp1_4f, _tmp2_4f );
                                _mm_store_ps(tmp_4f, _tmp_4f); 

                                _max_4f = _mm_max_ps(_score_4f[INST_Y], _tmp_4f);
                                _mm_store_ps(max_4f, _max_4f);

                                scoring_table_[iter_seq_x+0*t][iter_seq_y][INST_Y] = max_4f[0];
                                scoring_table_[iter_seq_x+1*t][iter_seq_y][INST_Y] = max_4f[1];
                                scoring_table_[iter_seq_x+2*t][iter_seq_y][INST_Y] = max_4f[2];
                                scoring_table_[iter_seq_x+3*t][iter_seq_y][INST_Y] = max_4f[3];
                                
//                                 scoring_table_[iter_seq_x+0*t][iter_seq_y][INST_Y] = std::max( {scoring_table_[iter_seq_x+0*t][iter_seq_y][INST_Y], tmp_4f[0] });
//                                 scoring_table_[iter_seq_x+1*t][iter_seq_y][INST_Y] = std::max( {scoring_table_[iter_seq_x+1*t][iter_seq_y][INST_Y], tmp_4f[1] });
//                                 scoring_table_[iter_seq_x+2*t][iter_seq_y][INST_Y] = std::max( {scoring_table_[iter_seq_x+2*t][iter_seq_y][INST_Y], tmp_4f[2] });
//                                 scoring_table_[iter_seq_x+3*t][iter_seq_y][INST_Y] = std::max( {scoring_table_[iter_seq_x+3*t][iter_seq_y][INST_Y], tmp_4f[3] });
                                
                                // update from up case
                                _tmp1_4f = _mm_sub_ps(_up_score_4f[MMATCH], _d_4f);
                                _tmp2_4f = _mm_sub_ps(_up_score_4f[INST_X], _e_4f);
                                _tmp_4f = _mm_max_ps( _tmp1_4f, _tmp2_4f );
                                _mm_store_ps(tmp_4f, _tmp_4f); 

                                _max_4f = _mm_max_ps(_score_4f[INST_X], _tmp_4f);
                                _mm_store_ps(max_4f, _max_4f);

                                scoring_table_[iter_seq_x+0*t][iter_seq_y][INST_X] = max_4f[0];
                                scoring_table_[iter_seq_x+1*t][iter_seq_y][INST_X] = max_4f[1];
                                scoring_table_[iter_seq_x+2*t][iter_seq_y][INST_X] = max_4f[2];
                                scoring_table_[iter_seq_x+3*t][iter_seq_y][INST_X] = max_4f[3];
                                
//                                 scoring_table_[iter_seq_x+0*t][iter_seq_y][INST_X] = std::max( {scoring_table_[iter_seq_x+0*t][iter_seq_y][INST_X], tmp_4f[0] });
//                                 scoring_table_[iter_seq_x+1*t][iter_seq_y][INST_X] = std::max( {scoring_table_[iter_seq_x+1*t][iter_seq_y][INST_X], tmp_4f[1] });
//                                 scoring_table_[iter_seq_x+2*t][iter_seq_y][INST_X] = std::max( {scoring_table_[iter_seq_x+2*t][iter_seq_y][INST_X], tmp_4f[2] });
//                                 scoring_table_[iter_seq_x+3*t][iter_seq_y][INST_X] = std::max( {scoring_table_[iter_seq_x+3*t][iter_seq_y][INST_X], tmp_4f[3] });

                        }

                        int levels_tobe_updated = 3;
                        for(int iter_interval = 1; iter_interval <= 3; iter_interval ++ ){
                                if(scoring_table_[iter_interval * t + 1][iter_seq_y] >= scoring_table_[iter_interval * t][iter_seq_y]){
                                        levels_tobe_updated -- ;
                                        break;
                                }
                        }

                        update_interval(levels_tobe_updated, scoring_table_, traceback_table_, iter_seq_y, t, _d_4f, _e_4f);

                        for(int iter_seq_x = 1; iter_seq_x <= seq_x.size(); iter_seq_x ++ ){
                                traceback_table_[iter_seq_x][iter_seq_y] = std::distance(
                                        scoring_table_[iter_seq_x][iter_seq_y].begin(),
                                        std::max_element(
                                                scoring_table_[iter_seq_x][iter_seq_y].begin(), 
                                                scoring_table_[iter_seq_x][iter_seq_y].end()
                                        )
                                );
                        }

                }

// TIMING END 
                auto __end__ = std::chrono::high_resolution_clock::now();
                auto __elapsed__ = std::chrono::duration_cast<std::chrono::nanoseconds>(__end__-__begin__);


                size_t tracer_x;
                size_t tracer_y;
                std::vector<float> largest3 (3, -1.0f);
                std::vector<std::vector<float>> largest3_idx (3, std::vector<float>(2,0));

                for(int i=1;i<seq_x_size_real;i++){
                        for(int j=1;j<seq_y.size();j++){
                                if( largest3[2] < scoring_table_[i][j][traceback_table_[i][j]] ){
                                        largest3[0] = largest3[1];
                                        largest3[1] = largest3[2];
                                        largest3[2] = scoring_table_[i][j][traceback_table_[i][j]];
                                        largest3_idx[0][0] = largest3_idx[1][0];
                                        largest3_idx[0][1] = largest3_idx[1][1];
                                        largest3_idx[1][0] = largest3_idx[2][0];
                                        largest3_idx[1][1] = largest3_idx[2][1];
                                        largest3_idx[2][0] = i;
                                        largest3_idx[2][1] = j;
                                }
                        }
                }
                                
                tracer_x = largest3_idx[2][0];
                tracer_y = largest3_idx[2][1];
                trace_and_align( seq_x, seq_y, tracer_x, tracer_y, scoring_table_, traceback_table_, seq_x_size_real);

                tracer_x = largest3_idx[1][0];
                tracer_y = largest3_idx[1][1];
                trace_and_align( seq_x, seq_y, tracer_x, tracer_y, scoring_table_, traceback_table_, seq_x_size_real);

                tracer_x = largest3_idx[0][0];
                tracer_y = largest3_idx[0][1];
                trace_and_align( seq_x, seq_y, tracer_x, tracer_y, scoring_table_, traceback_table_, seq_x_size_real);


                // output scoring table and trace back table
                std::ofstream scoring_file;
                scoring_file.open("./output/scoring_table");
                std::cout.precision(3);
                for(size_t iter_x=0; iter_x<seq_x.size()+1; iter_x++){
                        for(size_t iter_y=0; iter_y<seq_y.size()+1; iter_y++){
                                // std::cout << std::max({scoring_table_[iter_x][iter_y][MMATCH], scoring_table_[iter_x][iter_y][INST_X], scoring_table_[iter_x][iter_y][INST_Y]}) << "\t";
                                scoring_file << std::fixed << std::max({scoring_table_[iter_x][iter_y][MMATCH], scoring_table_[iter_x][iter_y][INST_X], scoring_table_[iter_x][iter_y][INST_Y]}) << "\t";
                        }
                        // std::cout << std::endl;
                        scoring_file << std::endl;
                }
                std::cout << "\n\n";

                std::ofstream traceback_file;
                scoring_file.open("./output/traceback_table");
                for(size_t iter_x=0; iter_x<seq_x.size()+1; iter_x++){
                        for(size_t iter_y=0; iter_y<seq_y.size()+1; iter_y++){
                                traceback_file << traceback_table_[iter_x][iter_y] << "\t";
                        }
                        traceback_file << std::endl;
                }


// OUTPUT TIMING RESULT
                std::cout << "Time Measured : \n";
                std::cout << __elapsed__.count() * 1e-9 << std::endl << std::endl;
                

// 		if (!hmm->isBasic()){
// 			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
// 			return;
// 		}
// 		
// 		//Initialize the traceback table
// 		if (traceback_table != NULL){
// 			delete traceback_table;
// 		}
// 		
// 		traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
// 		scoring_previous = new (std::nothrow) std::vector<float> (state_size,-INFINITY);
// 		scoring_current  = new (std::nothrow) std::vector<float> (state_size,-INFINITY);
// 		
// 		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
// 			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
// 			exit(2);
// 		}
// 		
// 		
// 		std::bitset<STATE_MAX> next_states;
// 		std::bitset<STATE_MAX> current_states;
// 		
// 		float  viterbi_temp(-INFINITY);
// 		float  emission(-INFINITY);
// 		bool	exDef_position(false);
// 		ending_viterbi_tb = -1;
// 		ending_viterbi_score = -INFINITY;
// 		
// 		state* init = hmm->getInitial();
// 		
// 		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
// 		std::bitset<STATE_MAX>* from_trans(NULL);
// 		
// 		//Calculate Viterbi from transitions from INIT (initial) state
// 		for(size_t st = 0; st < state_size; ++st){
// 			if ((*initial_to)[st]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
// 				
// 				viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) + getTransition(init, st, 0);
// 				
// 				if (viterbi_temp > -INFINITY){
// 					if ((*scoring_current)[st] < viterbi_temp){
// 						(*scoring_current)[st] = viterbi_temp;
// 					}
// 					next_states |= (*(*hmm)[st]->getTo());
// 				}
// 			}
// 		}
// 		
// 		//Each position in the sequence
// 		for(size_t position = 1; position < seq_size ; ++position ){
// 			
// 			//Swap current and previous viterbi scores
// 			scoring_previous->assign(state_size,-INFINITY);
// 			swap_ptr = scoring_previous;
// 			scoring_previous = scoring_current;
// 			scoring_current = swap_ptr;
// 			
// 			//Swap current_states and next states sets
// 			
// 			current_states.reset();
// 			current_states |= next_states;
// 			next_states.reset();
// 			
// 			if (exDef_defined){
// 				exDef_position = seqs->exDefDefined(position);
// 			}
// 			
// 			//Current states
// 			for (size_t st_current = 0; st_current < state_size; ++st_current){ //Current state that emits value
// 				
// 				//Check to see if current state is valid
// 				if (!current_states[st_current]){
// 					continue;
// 				}
// 				
// 				//Get emission of current state
// 				emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
// 				
// 				
// 				if (exDef_defined && exDef_position){
// 					emission += seqs->getWeight(position, st_current);
// 				}
// 
// 				
// 				if (emission == -INFINITY){
// 					continue;
// 				}
// 				
// 				//Get list of states that are valid previous states
// 				from_trans = (*hmm)[st_current]->getFrom();
// 				
// 				for (size_t st_previous = 0; st_previous < state_size ; ++st_previous) {  //for previous states
// 					if (!(*from_trans)[st_previous]){
// 						continue;
// 					}
// 					
// 					//Check that previous state has transition to current state
// 					//and that the previous viterbi score is not -INFINITY
// 					if ((*scoring_previous)[st_previous] != -INFINITY){
// 						viterbi_temp = getTransition((*hmm)[st_previous], st_current , position) + emission + (*scoring_previous)[st_previous];
// 						
// 						if (viterbi_temp > (*scoring_current)[st_current]){
// 							(*scoring_current)[st_current] = viterbi_temp;
// 							(*traceback_table)[position][st_current] = st_previous;
// 						}
// 						
// 						next_states |= (*(*hmm)[st_current]->getTo());
// 					}
// 				}
// 			}
// 		}
// 		
// 		//TODO:  Calculate ending and set the final viterbi and traceback pointer
// 		//Swap current and previous viterbi scores
// 		scoring_previous->assign(state_size,-INFINITY);
// 		swap_ptr = scoring_previous;
// 		scoring_previous = scoring_current;
// 		scoring_current = swap_ptr;
// 		
// 		
// 		//Calculate ending viterbi score and traceback from END state
// 		for(size_t st_previous = 0; st_previous < state_size ;++st_previous){
// 			if ((*scoring_previous)[st_previous] > -INFINITY){
// 				viterbi_temp = (*scoring_previous)[st_previous] + (*hmm)[st_previous]->getEndTrans();
// 				
// 				if (viterbi_temp > ending_viterbi_score){
// 					ending_viterbi_score = viterbi_temp;
// 					ending_viterbi_tb = st_previous;
// 				}
// 			}
// 		}
// 		
// 		delete scoring_previous;
// 		delete scoring_current;
// 		scoring_previous = NULL;
// 		scoring_current  = NULL;
	}
	
	
	void trellis::naive_viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		//TODO: determine which model and chose the type of algorithm to use;
		naive_viterbi();

	}
	
	
	void trellis::naive_viterbi(){
		traceback_table = new(std::nothrow) int_2D(seq_size, std::vector<int16_t>(state_size,-1));
		dbl_viterbi_score = new (std::nothrow) double_2D(seq_size, std::vector<double>(state_size,-INFINITY));
		
		double emission(-INFINITY);
		double viterbi_temp(-INFINITY);
		double trans(-INFINITY);
		double previous(-INFINITY);
		bool	exDef_position(false);
		ending_viterbi_tb = -1;
		ending_viterbi_score = -INFINITY;
		
		state* init = hmm->getInitial();
		
		//Calculate from Initial states
		for(size_t st = 0; st < state_size; ++st){
			viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) +  getTransition(init, st, 0);
			(*dbl_viterbi_score)[0][st]=viterbi_temp;
			(*traceback_table)[0][st]=-1;
		}
		
		//Calculate Forward for all states
		for (size_t position = 1 ; position < seq_size ; ++position){
			
			if (exDef_defined){
				exDef_position = seqs->exDefDefined(position);
			}
			
			for (size_t st_current = 0; st_current < state_size; ++st_current){
				//Calc emissions
				emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
				
				if (exDef_defined && exDef_position){
					emission += seqs->getWeight(position, st_current);
				}
				
				if (emission == -INFINITY){
					continue;
				}
				
				for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
					previous = (*dbl_viterbi_score)[position-1][st_previous];
					if (previous == -INFINITY){
						continue;
					}
					
					trans = getTransition(hmm->getState(st_previous), st_current, position);
					
					if (trans !=-INFINITY){
						viterbi_temp = emission + trans + previous;
						
						if (viterbi_temp > (*dbl_viterbi_score)[position][st_current]){
							(*dbl_viterbi_score)[position][st_current] = viterbi_temp;
							(*traceback_table)[position][st_current] = st_previous;
						}
					}
				}
			}
		}
		
		
		//Calculate Ending Transition
		ending_viterbi_tb = -1;
		for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
			if ((*hmm)[st_previous]->getEndTrans() != -INFINITY){
				if ((*dbl_viterbi_score)[seq_size-1][st_previous] != -INFINITY){
					viterbi_temp = (*dbl_viterbi_score)[seq_size-1][st_previous] + (*hmm)[st_previous]->getEndTrans();
					if (viterbi_temp > ending_viterbi_score){
						ending_viterbi_score = viterbi_temp;
						ending_viterbi_tb = st_previous;
					}
				}
			}
		}
		return;
	}
	
	
	
		
	
	
//	void trellis::viterbi(){
//		
//		//Initialize the traceback table
//		if (traceback_table != NULL){
//			delete traceback_table;
//		}
//		
//		traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
//		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
//		scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
//		
//		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
//			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
//			exit(2);
//		}
//		
//		
//		std::bitset<STATE_MAX> next_states;
//		std::bitset<STATE_MAX> current_states;
//		
//		double  viterbi_temp(-INFINITY);
//		double  emission(-INFINITY);
//		bool	exDef_position(false);
//		
//		
//		//If model is not a basic model, then we need to initialize the explicit duration vector
//		//bool extend_duration keeps track of whether the transition to same state was selected.
//		if (!hmm->isBasic()){
//			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
//			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
//		}
//		bool extend_duration(false);
//		std::vector<bool>* duration = hmm->get_explicit();
//		
//		state* init = hmm->getInitial();
//		
//		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
//		std::bitset<STATE_MAX>* from_trans(NULL);
//		
//		//Calculate Viterbi from transitions from INIT (initial) state
//		for(size_t i = 0; i < state_size; ++i){
//			if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
//				
//				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
//				
//				if (viterbi_temp > -INFINITY){
//					if ((*scoring_current)[i] < viterbi_temp){
//						(*scoring_current)[i] = viterbi_temp;
//					}
//					next_states |= (*(*hmm)[i]->getTo());
//				}
//			}
//		}
//		
//		//		for(size_t i=0; i < state_size; ++i){
//		//			std::cout << "Position: 0" << std::endl;
//		//			std::cout << exp((*viterbi_current)[i]) << std::endl;
//		//		}
//		//
//		
//		for(size_t position = 1; position < seq_size ; ++position ){
//			
//			//Swap current and previous viterbi scores
//			scoring_previous->assign(state_size,-INFINITY);
//			swap_ptr = scoring_previous;
//			scoring_previous = scoring_current;
//			scoring_current = swap_ptr;
//			
//			//Swap current_states and next states sets
//			
//			current_states.reset();
//			current_states |= next_states;
//			next_states.reset();
//			
//			if (exDef_defined){
//				exDef_position = seqs->exDefDefined(position);
//			}
//			
//			//TODO: Check use of external definitions below.
//			
//			std::cout << "\nPosition:\t" << position << "\n";
//			//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;
//			
//			for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
//				if (!current_states[i]){
//					continue;
//				}
//				
//				//current_state = (*hmm)[i];
//				//emission = current_state->get_emission(*seqs,position);
//				emission = (*hmm)[i]->get_emission_prob(*seqs, position);
//				
//				
//				//				std::cout << "State Emission:\t" << i << "\t" << exp(emission) << std::endl;
//				
//				if (exDef_defined && exDef_position){
//					emission += seqs->getWeight(position, i);
//				}
//				
//				from_trans = (*hmm)[i]->getFrom();
//				
//				for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
//					if (!(*from_trans)[j]){
//						continue;
//					}
//					
//					if ((*scoring_previous)[j] != -INFINITY){
//						viterbi_temp = getTransition((*hmm)[j], i , position) + emission + (*scoring_previous)[j];
//						
//						std::cout << exp(getTransition((*hmm)[j],i,position)) << std::endl;
//						
//						//						std::cout << "Temp Viterbi:\tTransFrom: "<< j << "\tto\t" << i << "\t" << viterbi_temp / log(2) << std::endl;
//						
//						
//						if (viterbi_temp > (*scoring_current)[i]){
//							//If transition is from same to same then if it is
//							//explicit duration we'll need to change
//							extend_duration = (i==j) ? true : false;
//							
//							(*scoring_current)[i] = viterbi_temp;
//							(*traceback_table)[position][i] = j;
//						}
//						
//						next_states |= (*(*hmm)[i]->getTo());
//					}
//				}
//				
//				//If explicit durations vector defined, and transition from-to same
//				//then we'll increment the value. Otherwise, set to zero;
//				if (explicit_duration_current){
//					if (extend_duration && (*duration)[i]){
//						(*explicit_duration_current)[i]=(*explicit_duration_previous)[i]+1;
//						extend_duration=false;
//					}
//					else{
//						(*explicit_duration_current)[i]=0;
//					}
//				}
//			}
//			
//			if(explicit_duration_current){
//				swap_ptr_duration = explicit_duration_previous;
//				explicit_duration_previous = explicit_duration_current;
//				explicit_duration_current= swap_ptr_duration;
//				explicit_duration_current->assign(state_size,0);
//			}
//			
//		}
//		
//		//TODO:  Calculate ending and set the final viterbi and traceback pointer
//		//Swap current and previous viterbi scores
//		scoring_previous->assign(state_size,-INFINITY);
//		swap_ptr = scoring_previous;
//		scoring_previous = scoring_current;
//		scoring_current = swap_ptr;
//		
//		for(size_t i = 0; i < state_size ;++i){
//			if ((*scoring_previous)[i] > -INFINITY){
//				viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
//				
//				if (viterbi_temp > ending_viterbi_score){
//					ending_viterbi_score = viterbi_temp;
//					ending_viterbi_tb = i;
//				}
//			}
//		}
//		
//		delete scoring_previous;
//		delete scoring_current;
//		scoring_previous = NULL;
//		scoring_current = NULL;
//	}
	
	
	//TODO:  Need to test and finalize these algorithms perform
	/*	Need to simplify the basic calls so that it checks model and chooses the
	 algorithm to perform
	 
	 Need to establish duration algorithms for forward/backward that first to viterbi
	 to calculate the transition probability.
	 
	 */
	
	
	//! Sparse Complex Viterbi
	//! Stores the transition duration probabilities in a hashmap (memory efficient, slower)
	//! The duratio probabilities can then be used in forward and backward algorithms.
	void trellis::sparse_complex_viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,1));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
		double  transition_prob(-INFINITY);
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = -1;
		
		
		//If model is not a basic model, then we need to initialize the explicit duration vector
		//bool extend_duration keeps track of whether the transition to same state was selected.
		if (!hmm->isBasic()){
			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
		}
		
		bool extend_duration(false);
		
		// Get list of States with explicit duration
		std::vector<bool>* duration = hmm->get_explicit();
		
		//Initialize Duration storage table
		complex_transitions = new std::vector<std::vector< std::map<uint16_t,double>* >* > (state_size,NULL);
		for(size_t i=0; i < state_size; i++){
			if((*duration)[i]){
				(*complex_transitions)[i] = new std::vector<std::map<uint16_t,double>* > (seq_size, NULL);
			}
		}
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				//Transitions here are guarenteed to be standard from the initial state
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (viterbi_temp > -INFINITY){
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states |= (*(*hmm)[i]->getTo());
                }
            }
        }
		
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
            //Swap current and previous viterbi scores
            scoring_previous->assign(state_size,-INFINITY);
            swap_ptr = scoring_previous;
			scoring_previous = scoring_current;
			scoring_current = swap_ptr;
            
            //Swap current_states and next states sets
			
			current_states.reset();
            current_states |= next_states;
            next_states.reset();
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
			
			//TODO: Check use of external definitions below.
			
			//std::cout << "\nPosition:\t" << position << "\n";
			//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;
			
            for (size_t st_current = 0; st_current < state_size; ++st_current){ //i is current state that emits value
                if (!current_states[st_current]){
                    continue;
                }
				
                emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
				
				
				//Check External definitions
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, st_current);
                }
                
				//Get list of state that transition to current state
				from_trans = (*hmm)[st_current]->getFrom();
				
                for (size_t st_previous = 0; st_previous < state_size ; ++st_previous) {  //j is previous state
                    if (!(*from_trans)[st_previous]){ //if transition not possible next
                        continue;
                    }
					
                    if ((*scoring_previous)[st_previous] != -INFINITY){
						
						transition_prob = getTransition((*hmm)[st_previous], st_current , position);
						
						if ((*duration)[st_previous]){
							if ((*(*complex_transitions)[st_current])[position] == NULL){
								(*(*complex_transitions)[st_current])[position] = new std::map<uint16_t,double>;
							}
							
							(*(*(*complex_transitions)[st_current])[position])[st_previous] = transition_prob;
						}
						
                        viterbi_temp = transition_prob + emission + (*scoring_previous)[st_previous];
                        
						
						if (viterbi_temp > (*scoring_current)[st_current]){
							//If transition is from same to same then if it is
							//explicit duration we'll need to change
							extend_duration = (st_current==st_previous) ? true : false;
							
                            (*scoring_current)[st_current] = viterbi_temp;
                            (*traceback_table)[position][st_current] = st_previous;
                        }
						
						next_states |= (*(*hmm)[st_current]->getTo());
                    }
                }
				
				//If explicit durations vector defined, and transition from-to same
				//then we'll increment the value. Otherwise, set to zero;
				if (explicit_duration_current){
					if (extend_duration && (*duration)[st_current]){
						if ((*explicit_duration_previous)[st_current]==0){
							(*explicit_duration_current)[st_current]=2;
						}
						else{
							(*explicit_duration_current)[st_current]=(*explicit_duration_previous)[st_current]+1;
						}
						extend_duration=false;
					}
					else{
						(*explicit_duration_current)[st_current]=0;
					}
				}
            }
			
			if(explicit_duration_current){
				swap_ptr_duration = explicit_duration_previous;
				explicit_duration_previous = explicit_duration_current;
				explicit_duration_current= swap_ptr_duration;
				explicit_duration_current->assign(state_size,0);
			}
            
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
        
        for(size_t st_previous = 0; st_previous < state_size ;++st_previous){
            if ((*scoring_previous)[st_previous] > -INFINITY){
                viterbi_temp = (*scoring_previous)[st_previous] + (*hmm)[st_previous]->getEndTrans();
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = st_previous;
                }
            }
        }
        
        delete scoring_previous;
        delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
	
	//!Store transition in a table  (more memory, but faster than sparse complex)
	//!Need testing and more develepment;
	void trellis::fast_complex_viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
		ending_viterbi_tb = -1;
		ending_viterbi_score = -INFINITY;
		
		
		//If model is not a basic model, then we need to initialize the explicit duration vector
		//bool extend_duration keeps track of whether the transition to same state was selected.
		if (!hmm->isBasic()){
			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
		}
		bool extend_duration(false);
		std::vector<bool>* duration = hmm->get_explicit();
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (viterbi_temp > -INFINITY){
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states |= (*(*hmm)[i]->getTo());
//					if ((*duration)[i]){
//						(*explicit_duration_previous)[i]=1;
//					}
                }
            }
        }
		
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
            //Swap current and previous viterbi scores
            scoring_previous->assign(state_size,-INFINITY);
            swap_ptr = scoring_previous;
			scoring_previous = scoring_current;
			scoring_current = swap_ptr;
            
            //Swap current_states and next states sets
			
			current_states.reset();
            current_states |= next_states;
            next_states.reset();
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
			
			//TODO: Check use of external definitions below.
			
//			std::cout << "\nPosition:\t" << position << "\n";
//			std::cout << "Letter:\t" << seqs->seqValue(0, position)+1 << std::endl;
			
            for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
                if (!current_states[i]){
                    continue;
                }
				
                //current_state = (*hmm)[i];
                //emission = current_state->get_emission(*seqs,position);
                emission = (*hmm)[i]->get_emission_prob(*seqs, position);
				
				
//				std::cout << "State Emission:\t" << i << "\t" << exp(emission) << std::endl;
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
				from_trans = (*hmm)[i]->getFrom();
				
                for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
                    if (!(*from_trans)[j]){
                        continue;
                    }
					
                    if ((*scoring_previous)[j] != -INFINITY){
                        viterbi_temp = getTransition((*hmm)[j], i , position) + emission + (*scoring_previous)[j];
                        
//						std::cout << "TransFrom: "<< j << "\tto\t" << i << "\t" << exp(getTransition((*hmm)[j],i,position)) << std::endl;
//						
//						std::cout << "Temp Viterbi:\t" << viterbi_temp << std::endl;
                        
						
						if (viterbi_temp > (*scoring_current)[i]){
							//If transition is from same to same then if it is
							//explicit duration we'll need to change
							extend_duration = (i==j) ? true : false;
							
                            (*scoring_current)[i] = viterbi_temp;
                            (*traceback_table)[position][i] = j;
                        }
						
						next_states |= (*(*hmm)[i]->getTo());
                    }
                }
				
				//If explicit durations vector defined, and transition from-to same
				//then we'll increment the value. Otherwise, set to zero;
				if (explicit_duration_current){
					if (extend_duration && (*duration)[i]){
						if ((*explicit_duration_previous)[i]==0){
							(*explicit_duration_current)[i]=2;
						}
						else{
							(*explicit_duration_current)[i]=(*explicit_duration_previous)[i]+1;
						}
						extend_duration=false;
					}
					else{
						(*explicit_duration_current)[i]=0;
					}
				}
            }
			
			if(explicit_duration_current){
				swap_ptr_duration = explicit_duration_previous;
				explicit_duration_previous = explicit_duration_current;
				explicit_duration_current= swap_ptr_duration;
				explicit_duration_current->assign(state_size,0);
			}
            
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] > -INFINITY){
                viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                }
            }
        }
        
        delete scoring_previous;
        delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
	

}

