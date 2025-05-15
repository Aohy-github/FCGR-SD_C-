
#include "readSeq.hpp"

std::vector<Seq> get_sequence(const std::string& filename , int K){
    int max_lens = 0, min_lens = INT32_MAX , mid_lens = 0;
    std::ifstream infile(filename);
    std::vector<Seq> seqs;
    if(!infile.is_open()) {
        printf("open file failed!");
        exit(0);
    }
    else {
        std::string line;
        std::string seqName;
        std::string content;
        while(getline(infile, line)){
            if(line[0] == '>') {
                if(content.size() > 0) {

                    Seq seq(seqName, content, K);
                    seqs.push_back(seq);
                }
                seqName = line;
                content.clear();
            }else {
                content += line;
            }
        }
        if(content.size() > 0) {
            Seq seq(seqName, content, K);

            seqs.push_back(seq);
        }

        infile.close();
    }
    return seqs;
}
