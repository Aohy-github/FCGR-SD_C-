#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <boost/program_options.hpp>
#include <sys/stat.h>
#include <iomanip>
namespace po = boost::program_options;

#include "readSeq.hpp"
#include "Seq.hpp"
#include "fcgr.hpp"
#include "ThreadPool.h"

int main(int argc , char *argv[])
{
    std::string filename;
    int kmerSize = 6;
    int N_thread = 1;
    double SVD_P = 0.0;
    po::options_description desc("选项:");
    desc.add_options()

    ("help,h", "显示帮助信息")
    ("input,i", po::value<std::string>(&filename)->required(), "输入 FASTA 文件路径")
    ("kmerSizem,k", po::value<int>(&kmerSize)->default_value(6), "kmer大小")
    ("SVD_P,p", po::value<double>(&SVD_P)->default_value(0.75), "kmer大小")
    ("thread,t", po::value<int>(&N_thread)->default_value(1), "线程数量");

    po::variables_map vm;
    try {
        // 解析命令行参数
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }

        // 检查必需参数
        po::notify(vm);
    }
    catch(po::error& e) {
        std::cerr << "error_massage: " << e.what() << "\n";
        std::cerr << desc << "\n";
        return 1;
    }

    std::cout << "输入文件: " << filename << "\n";
    // 读取序列
    std::vector<Seq> seqs = get_sequence(filename, kmerSize);
    FCGR fcgr;
    std::cout << "N_thread : " << N_thread << std::endl;
    ThreadPool threadPool(N_thread);
    
    std::vector<std::future<Eigen::MatrixXd>> res_fcgr_list;
    std::vector<std::future<Eigen::MatrixXd>> res_SVD_list;

    struct Com_fcgr_task{
        FCGR fcgr;
        std::string content;
        int kmerSize;

        Eigen::MatrixXd operator()() const{
            return fcgr.computeMatrix(content, kmerSize);
        }
    };
    struct Com_SVD_task{
        FCGR fcgr;
        Eigen::MatrixXd matrix;
        double P;

        Eigen::MatrixXd  operator()() const{
            return fcgr.computerSVD(matrix , P);
        }
    };
    // 计算FCGR矩阵
    for (auto& seq : seqs) {
        Com_fcgr_task task{fcgr , seq.getSeqContent() , kmerSize};
        std::future<Eigen::MatrixXd> res_fcgr = threadPool.enqueue(task);
        
        res_fcgr_list.emplace_back(std::move(res_fcgr));
        
    }
    if (res_fcgr_list.size() != seqs.size()) {
        std::cerr << "Mismatch in task enqueuing, program may hang!" << std::endl;
        exit(0);
    }

    for (size_t i = 0; i < res_fcgr_list.size(); ++i) {
        seqs[i].setMatrix(res_fcgr_list[i].get());
    }
    
    int count_temp = 0;
    std::cout << "seqs.size() :" << seqs.size() << std::endl;
    for(auto& one: seqs){
        Eigen::MatrixXd matrix = one.getMatrix();
        Com_SVD_task task{fcgr , matrix , SVD_P};
        auto res_SVD = threadPool.enqueue(task);
        res_SVD_list.emplace_back(std::move(res_SVD));
    }

    std::cout << "计算完成\n";
    for(size_t i = 0; i < res_SVD_list.size() ; i++){
        seqs[i].setrankMatrix(res_SVD_list[i].get());
    }

    
    // 计算Grassmann距离
    int n = seqs.size();
    //std::cout << "序列数量: " << n << "\n";
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n, n);
    for(int i = 0 ; i < n; i++){
        for(int j = i + 1 ; j < n; j++){
            Eigen::MatrixXd rankMatrix1 = seqs[i].getRankMatrix();
            Eigen::MatrixXd rankMatrix2 = seqs[j].getRankMatrix();
            D(i,j) = fcgr.computeGrassmannDistance(rankMatrix1, rankMatrix2);
            D(j,i) = D(i,j);
        }
    }
    std::string dis_filename = "distance_matrix.txt";
    std::ofstream outfile(dis_filename);
    if(!outfile.is_open()) {
        printf("open file failed");
        exit(0);
    }
    outfile << seqs.size() << std::endl;
    for(int i = 0 ; i < seqs.size() ; i ++){
        outfile <<std::left << std::setw(15) << seqs[i].getSeqName().substr(1);
        for(int j = 0; j < seqs.size() ; j ++){
            outfile <<std::left << std::setw(10) << D(i,j);
        }
        outfile << std::endl;
    }

    outfile.close();



     std::string cmnd = "./decenttree ";
        std::string out_res_file = "out.tree";
        cmnd.append("-in ").append("distance_matrix.txt ")
            .append("-t BIONJ ")
            .append("-nt 16 ")
            .append("-no-banner ")
            .append("-out ").append(out_res_file).append("> /dev/null");

        int res = system(cmnd.c_str());
        if (res != 0){
            std::cout << "run align softward false";
            exit(0);
        }else{
            std::cout << "Tree file is ok: " << out_res_file;
        }
    return 0;
}
