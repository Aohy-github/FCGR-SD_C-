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


int main(int argc , char *argv[])
{
    std::string filename;
    int kmerSize = 6;
    po::options_description desc("选项:");
    desc.add_options()

    ("help,h", "显示帮助信息")
    ("input,i", po::value<std::string>(&filename)->required(), "输入 FASTA 文件路径")
    ("kmerSizem,k", po::value<int>(&kmerSize)->default_value(6), "kmer大小");

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
    std::vector<Seq> seqs = get_sequence(filename, kmerSize);
    FCGR fcgr;
    // for(int i =0 ; i < seqs.size(); i++) {
    //     std::cout << "序列名称: " << seqs[i].getSeqName() << "\n";
    //     std::cout << "序列内容: " << seqs[i].getSeqContent().size() << "\n";
    // }

    // 计算FCGR矩阵
    for (auto& seq : seqs) {
        seq.setMatrix(fcgr.computeMatrix(seq.getSeqContent() , kmerSize));
        
        Eigen::MatrixXd matrix = seq.getMatrix();
        
        std::cout << "序列名称: " << seq.getSeqName() << "\n";
        std::cout << "FCGR矩阵:\n" << matrix << "\n";
        
        Eigen::MatrixXd rankMatrix = seq.getRankMatrix();
        seq.setrankMatrix(fcgr.computerSVD(seq.getMatrix()));
        std::cout << "秩矩阵:\n" << rankMatrix << "\n";
    }
    std::cout << "计算完成\n";
    // 计算Grassmann距离
    int n = seqs.size();
    std::cout << "序列数量: " << n << "\n";
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
