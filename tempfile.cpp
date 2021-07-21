//Each position in the sequence
for(size_t iter_seq_x = 1; iter_seq_x < seq_x_size+1; ++iter_seq_x){
        for(size_t iter_seq_y = 1; iter_seq_y < seq_y_size+1; ++iter_seq_y){

                scoring_table[iter_seq_x][iter_seq_y][MMATCH] = 
                        hmm->alignLogProb[seqs_x[iter_seq_x]][seqs_y[iter_seq_y]] + std::max({
                                (1- 2*hmm->delta - hmm->tau) * scoring_table[iter_seq_x-1][iter_seq_y-1][MMATCH],
                                (1- hmm->epsilon - hmm->tau) * scoring_table[iter_seq_x-1][iter_seq_y-1][INST_X],
                                (1- hmm->epsilon - hmm->tau) * scoring_table[iter_seq_x-1][iter_seq_y-1][INST_Y]});
                                
                scoring_table[iter_seq_x][iter_seq_y][INST_X] =
                        std::max({
                                hmm->delta * scoring_table[iter_seq_x-1][iter_seq_y][MMATCH] - hmm->d,
                                hmm->epsilon * scoring_table[iter_seq_x-1][iter_seq_y][INST_X] - hmm->e});
                                
                scoring_table[iter_seq_x][iter_seq_y][INST_Y] =
                        std::max({
                                hmm->delta * scoring_table[iter_seq_x][iter_seq_y-1][MMATCH] - hmm->d,
                                hmm->epsilon * scoring_table[iter_seq_x][iter_seq_y-1][INST_Y]} - hmm->e);

                traceback_table[iter_seq_x-1][iter_seq_y-1] = std::distance( 
                        std::max_element( scoring_table[iter_seq_x-1][iter_seq_y-1].begin(), 
                        std::max_element( scoring_table[iter_seq_x-1][iter_seq_y-1].end())));
        }
}

// Termination
final_score = std::max({
        scoring_table[seq_x_size][seq_y_size][MMATCH],
        scoring_table[seq_x_size][seq_y_size][INST_X] + hmm->c,
        scoring_table[seq_x_size][seq_y_size][INST_Y] + hmm->c});
