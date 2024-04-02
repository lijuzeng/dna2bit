//
// Created by ljz on 24-1-2.
//

#ifndef DNA_TO_BIT_COMMAND_HPP
#define DNA_TO_BIT_COMMAND_HPP

#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
using namespace std;

void MainUsage(const char *const *argv) {
    cerr << "Usage:\n\t" << argv[0] << " <command> [options] [arguments ...]\n\n"
         << "command:\n\n"
         << "  sketch\tCreate sketch of input sequence file\n\n"
         << "  dist\t\tEstimate distance between sketched files\n\n"
         << "Note:\n\n"
         << "  To see command-specific usage, please enter\n"
         << "    " << argv[0] << " <command> --help\n";
    exit(0);
}

class SketchArgs{
public:
    vector<string> file_names_list;
    string list_files_name;
    size_t kmer_len = 17;
    size_t bit_len = 8192;
    size_t nthreads = 0;
    size_t hash_type = 0;
    SketchArgs() = default;
    ~SketchArgs() = default;
    void Usage(const char *const *argv);
    void GetArgs(const int argc, char *const *argv);
};

void SketchArgs::Usage(const char *const *argv) {
    cerr << "Usage : " << argv[0] << " " << argv[1] << " [options] [arguments ...]\n\n";
    cerr << "Arguments:\n\n"
         << "  One or more filenames. If zero filename, then read from each line in fnameslist. \n"
         << "  Each line of fnameslist is a path to a sequence file\n\n";
    cerr << "Option with [default values] : \n\n";
    cerr << "  --help / -h : Show this help message.\n\n";
    cerr << "  --kmer_len / -k : K-mer length used to generate hash value. [" << kmer_len << "]\n\n";
    cerr << "  --bit_len / -l : The length of generated bit of sequence file. [" << bit_len << "]\n\n";
    cerr << "  --nthreads / -n : This many threads will be spawned for processing. \n"
            "                    if this option be chosen, the option F must be chosen[" << nthreads << "]\n\n";
    cerr << "  --listfname / -F : In this file each line have the paths of genome files. [" << list_files_name << "]\n\n";
    cerr << "  --hash_type / -t : Type of hash function. \n"
         << "   0 : wyhash. \n"
         << "   1 : rolling hash. \n"
         << "   2 : murmur hash. [" << hash_type << "]\n\n";
    cerr << "Note : \n\n"
         << "  For general usage, please enter \n"
         << "    " << argv[0] << " --help\n\n";
    exit(0);
}

void SketchArgs::GetArgs(const int argc, char *const *argv) {
    const char *short_option = "hk:l:n:F:t:？";
    const option long_option[] = {
            {"help", no_argument, nullptr, 'h'},
            {"kem_len", required_argument, nullptr, 'k'},
            {"bit_len", required_argument, nullptr, 'l'},
            {"nthreads", required_argument, nullptr, 'n'},
            {"listfname", no_argument, nullptr, 'F'},
            {"hash_type", required_argument, nullptr, 't'}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, short_option, long_option, nullptr)) != -1) {
        switch (opt) {
            case 'h' : Usage(argv); break;
            case 'k' : kmer_len = atoi(optarg);     break;
            case 'l' : bit_len = atoi(optarg);      break;
            case 'n' : nthreads = atoi(optarg);     break;
            case 'F' : list_files_name = optarg;          break;
            case 't' : hash_type = atoi(optarg);    break;
            case '?' : default : Usage(argv); break;
        }
    }
    if (list_files_name.empty() && argc != optind + 2) {
        cerr << "incorrect usage of sketch with one file. \n"
             << "Enter \"--help\" (without quotes) to show usage. " << endl;
        for (string line; getline(cin, line);) {
            if (line == "--help") { Usage(argv); }
        }
    }
    if (!list_files_name.empty() && argc != optind + 1) {
        cerr << "incorrect usage of sketch with files list. \n"
             << "Enter \"--help\" (without quotes) to show usage. " << endl;
        for (string line; getline(cin, line);) {
            if (line == "--help") { Usage(argv); }
        }
    }
    if (!list_files_name.empty()) {
        char tmp[1024];
        FILE *p = fopen(list_files_name.c_str(), "r");
        if (p == nullptr) {
            printf("%s unable open\n", list_files_name.c_str());
            exit(0);
        }
        while (fgets(tmp, 1024, p) != nullptr) {file_names_list.push_back(tmp);}
        size_t files_size = file_names_list.size();
        for (size_t i = 0; i < files_size; ++i) { file_names_list[i].pop_back(); }
        fclose(p);
    }
}

void GetKandL(const char *str, size_t &kmer, size_t &bit) {
    char *temp_k, *temp_l;
    size_t len = strlen(str), loc0 = 0, loc1 = 0, loc2 = 0, loc3 = 0;
    for (size_t i = 0; i < len; ++i) {
        if (strncmp(str + i, ".k.", 3) == 0) {loc0 = i + 3;}
        if (strncmp(str + i, ".l.", 3) == 0) {loc1 = i; loc2 = i + 3;}
        if (strncmp(str + i, ".bit", 4) == 0) {loc3 = i;}
    }
    if (loc0 == 0 && loc1 == 0 && loc2 == 0 && loc3 == 0) {
        cerr << "The input file is incorrect." << endl;
        exit(0);
    }
    temp_k = (char *)malloc((loc1 - loc0) * sizeof(char));
    temp_l = (char *)malloc((loc3 - loc2) * sizeof(char));
    strncpy(temp_k, str + loc0, loc1 - loc0);
    strncpy(temp_l, str + loc2, loc3 - loc2);
    kmer = atoi(temp_k);
    bit = atoi(temp_l);
    free(temp_k);
    free(temp_l);
}

class DistArgs{
public:
    vector<string> file_names_list;
    string list_files_name;
    string output_filename = "distance.csv";
    size_t nthreads = 0;
    size_t kmer_len, bit_len;
    bool triangular = true;
    DistArgs() = default;
    ~DistArgs() = default;
    void Usage(const char *const *argv);
    void GetArgs(const int argc, char *const *argv);
};

void DistArgs::Usage(const char *const *argv) {
    cerr << "Usage : " << argv[0] << ' ' << argv[1] << " [options] [argument ... ]\n\n";
    cerr << "Argument : \n\n"
         << "  One bit files list or two bit files, if one list, then read each line in the list. \n"
         << "  Each line of the list is a bit file generated from sketch command. \n\n";
    cerr << "option with [default value] : \n\n";
    cerr << "  --help / -h : Show this help message. \n\n";
    cerr << "  --pair / -p : Input two bit file, and calculate distance metric(xor->popcount). \n\n";
    cerr << "  --listfname / -F : In this file each line have the paths of bit files. [" << list_files_name << "]\n\n";
    cerr << "  --nthreads / -n : This many threads will be spawned for processing. \n"
            "                    if this option be chosen, the option F must be chosen[" << nthreads << "]\n\n";
    cerr << "  --whole_matrix / -w : The default output martix is upper triangular.\n "
            "                        If this option be chosen the output is whole matrix.\n\n";
    cerr << "  --outfname / -o : Generate a upper triangular matrix with distance between each bit files. \n"
         << "                    if this option be chosen, the option F must be chosen. [" << output_filename << "]\n\n";
    cerr << "Note : \n\n";
    cerr << "  For general usage, please enter \n";
    cerr << "    " << argv[0] << " --help / -h \n\n";
    exit(0);
}

void DistArgs::GetArgs(const int argc, char *const *argv) {
    const char *short_option = "pF:n:wo:?";
    const option long_option[] = {
            {"help", no_argument, nullptr, 'h'},
            {"pair", no_argument, nullptr, 'p'},
            {"listfname", required_argument, nullptr, 'F'},
            {"nthreads", required_argument, nullptr, 'n'},
            {"whole_matrix", no_argument, nullptr, 'w'},
            {"outfname", required_argument, nullptr, 'o'}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, short_option, long_option, nullptr)) != -1) {
        switch (opt) {
            case 'h' : Usage(argv);                 break;
            case 'F' : list_files_name = optarg;    break;
            case 'p' : break;
            case 'n' : nthreads = atoi(optarg); break;
            case 'w' : triangular = false;          break;
            case 'o' : output_filename = optarg;    break;
            case '?' : default : Usage(argv);       break;
        }
    }
    if (list_files_name.empty() && argc != optind + 3) {
        cerr << "incorrect usage of dist with two bit files. \n"
             << "Enter \"--help\" (without quotes) to show usage. " << endl;
        for (string line; getline(cin, line);) {
            if (line == "--help") { Usage(argv); }
        }
    }
    if (!list_files_name.empty() && argc != optind + 1) {
        cerr << "incorrect usage of dist with bit files list. \n"
             << "Enter \"--help\" (without quotes) to show usage. " << endl;
        for (string line; getline(cin, line);) {
            if (line == "--help") { Usage(argv); }
        }
    }
    if (!list_files_name.empty()) {
        char tmp[1024];
        FILE *p = fopen(list_files_name.c_str(), "r");
        if (p == nullptr) {
            printf("%s unable open\n", list_files_name.c_str());
            exit(0);
        }
        while (fgets(tmp, 1024, p) != nullptr) {file_names_list.push_back(tmp);}
        size_t files_size = file_names_list.size();
        for (size_t i = 0; i < files_size; ++i) { file_names_list[i].pop_back();}
        GetKandL(file_names_list[0].c_str(), kmer_len, bit_len);
        fclose(p);
    } else {
        GetKandL(argv[optind + 1], kmer_len, bit_len);
    }
}

#endif //DNA_TO_BIT_COMMAND_HPP
