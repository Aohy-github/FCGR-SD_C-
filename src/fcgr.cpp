#include "fcgr.hpp"
#include <cmath>
#include <random>
#include <iostream>
#include "spdlog/spdlog.h"
FCGR::FCGR() {
    std::cout << "FCGR constructor called" << std::endl;
}
Eigen::MatrixXd FCGR::computeMatrix(const std::string& seq ,const int& k) {
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(1 << k, 1 << k);
    int len = seq.length();
    if (len < k){
        std::cout << "序列长度小于k，返回零矩阵" << std::endl;
        exit(0);
        return matrix; // 如果序列长度小于k，返回零矩阵
    } 

    for (int i = 0; i <= len - k; ++i) {
        std::string kmer = seq.substr(i, k);
        auto [x, y] = kmerToIndex(kmer , k);
        if (x >= 0 && y >= 0)
            matrix(x, y) += 1;
    }
    return matrix;
}
char otherCharToATCG(char num) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, 3); // 用于随机选择 A、T、C、G

    // 用于从不同碱基集合中随机选择
    int arr_ATCG[] = {'A', 'T', 'C', 'G'};
    int arr_atcg[] = {'a', 't', 'c', 'g'};
    int arr_AG[] = {'A', 'G'};
    int arr_CT[] = {'C', 'T'};
    int arr_CG[] = {'C', 'G'};
    int arr_AT[] = {'A', 'T'};
    int arr_GT[] = {'G', 'T'};
    int arr_AC[] = {'A', 'C'};
    int arr_CGT[] = {'C', 'G', 'T'};
    int arr_AGT[] = {'A', 'G', 'T'};
    int arr_ACT[] = {'A', 'C', 'T'};
    int arr_ACG[] = {'A', 'C', 'G'};
    
    switch(num) {
        case 'N': return arr_ATCG[dis(gen)];
        case 'R': return arr_AG[dis(gen) % 2]; // A or G
        case 'Y': return arr_CT[dis(gen) % 2]; // C or T
        case 'S': return arr_CG[dis(gen) % 2]; // G or C
        case 'W': return arr_AT[dis(gen) % 2]; // A or T
        case 'K': return arr_GT[dis(gen) % 2]; // G or T
        case 'M': return arr_AC[dis(gen) % 2]; // A or C
        case 'B': return arr_CGT[dis(gen) % 3]; // C or G or T
        case 'D': return arr_AGT[dis(gen) % 3]; // A or G or T
        case 'H': return arr_ACT[dis(gen) % 3]; // A or C or T
        case 'V': return arr_ACG[dis(gen) % 3]; // A or C or G
        case 'U': return 'T'; // U 转化为 T（用于 RNA 中）
        case 'u': return 't'; // U 转化为 T（用于 RNA 中）
        // case '-': return '-'; // gap
        case ':': return arr_ATCG[dis(gen)];
        // 处理有效的 ATCG 字符
        case 'A':
        case 'T':
        case 'C':
        case 'G':
        case 'a':
        case 't':
        case 'c':
        case 'g': return num;

        // 默认处理其他字符
        default: {
            std::cerr << "Unknown character: " << char(num) << std::endl;
            return arr_ATCG[dis(gen)];
            // exit(EXIT_FAILURE);
        }
    }
}

int FCGR::baseToInt(char base) {
    switch (otherCharToATCG(base)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'a': return 0;
        case 't': return 1;
        case 'c': return 2;
        case 'g': return 3;
        default: return -1;
    }
}

std::pair<int, int> FCGR::kmerToIndex(const std::string& kmer , const int& k) {
    int x = 0, y = 0;
    for (int i = 0; i < k; ++i) {
        int b = baseToInt(kmer[i]);
        if (b == -1) return {-1, -1}; 
        int bit = 1 << (k - i - 1);
        x |= ((b >> 1) & 1) * bit;
        y |= (b & 1) * bit;
    }
    return {x, y};
}


Eigen::MatrixXd FCGR::computerSVD(Eigen::MatrixXd Matrix, double P){
    //spdlog::info("=== start new ComSVD === ");
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    const auto& S = svd.singularValues();   // 奇异值 σ1, σ2, ..., σn
    const auto& U = svd.matrixU();          // 左奇异向量

    // 计算总能量
    double total = S.sum();

    // 寻找前 r 个奇异值使其累加占比达到 P
    double acc = 0.0;
    int r = 0;
    for (; r < S.size(); ++r) {
        acc += S[r];
        if (acc / total >= P) break;
    }
    r += 1;  // 因为索引从0开始，实际需要 r+1 列

    if (r > U.cols()) {
        std::cerr << "Error: r exceeds the number of columns in U." << std::endl;
        exit(EXIT_FAILURE);
        // return Eigen::MatrixXd(); // 错误处理
    }
    //spdlog::info("=== End new ComSVD === ");
    return U.leftCols(r);
}

double FCGR::computeGrassmannDistance(const Eigen::MatrixXd A, const Eigen::MatrixXd B, double threshold) {
    if (A.rows() != B.rows()) {
        std::cerr << "Error: The number of rows in A and B must be the same." << std::endl;
        return -1; // 错误处理
    }
    Eigen::MatrixXd M = A.transpose() * B;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd sigma = svd.singularValues();

    int r = sigma.size();
    Eigen::VectorXd theta = Eigen::VectorXd::Zero(r);  // 安全初始化

    for (int i = 0; i < r; ++i) {
        double cos_theta = std::clamp(sigma(i), -1.0, 1.0);

        if (cos_theta < threshold) {
            theta(i) = std::acos(cos_theta);
        } else {
            Eigen::VectorXd v = svd.matrixV().col(i);              // 确保是向量类型
            Eigen::VectorXd v_cos = B * v;
            Eigen::VectorXd proj = A * (A.transpose() * v_cos);
            double sin_val = (v_cos - proj).norm();
            sin_val = std::clamp(sin_val, 0.0, 1.0);
            theta(i) = std::asin(sin_val);
        }
    }

    return std::sqrt(theta.squaredNorm());
}