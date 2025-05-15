#ifndef SEQ_HPP
#define SEQ_HPP

#include <string>
#include "../include/fcgr.hpp"
class Seq{
    public :
        Seq() = default;
        Seq(const std::string& seqName, const std::string& content , int kmerLength): 
                    seqName(seqName) , 
                    SeqContent(content),
                    k(kmerLength) {}

        ~Seq() = default;
        std::string getSeqName() const { return seqName; }
        std::string getSeqContent() const { return SeqContent; }
        Eigen::MatrixXd getMatrix() {
            return matrix;
        }
        Eigen::MatrixXd getRankMatrix() const {
            return rankMatrix;
        }
        void setrankMatrix(const Eigen::MatrixXd rankMatrix) {
            this->rankMatrix = rankMatrix;
        }
        void setMatrix(const Eigen::MatrixXd matrix) {
            this->matrix = matrix;
        }

    private:
        std::string seqName;
        std::string SeqContent;
        Eigen::MatrixXd matrix;
        Eigen::MatrixXd rankMatrix;
        int k;
};


#endif