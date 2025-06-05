#ifndef FCGR_HPP
#define FCGR_HPP

#include <string>
#include <unordered_map>
#include <Eigen/Dense>
#include <iostream>
class FCGR{
    public :
        FCGR();
        ~FCGR() = default;
        static Eigen::MatrixXd computeMatrix(const std::string& seq ,const int& k);
        
        static Eigen::MatrixXd computerSVD(Eigen::MatrixXd Matrix, double P = 0.9);
        static double computeGrassmannDistance(const Eigen::MatrixXd A, const Eigen::MatrixXd B, double threshold = std::cos(1e-5));
        
        static std::pair<int, int> kmerToIndex(const std::string& kmer , const int& k);
        static int baseToInt(char base);
};
#endif