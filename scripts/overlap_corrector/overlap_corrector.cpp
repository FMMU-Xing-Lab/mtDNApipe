#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

// 基础数据结构
class Read {
public:
    string* mSeq;
    string* mQuality;
    string* mName;

    Read(string* name, string* seq, string* quality) {
        mName = name;
        mSeq = seq;
        mQuality = quality;
    }

    int length() { return mSeq->length(); }
};

class OverlapResult {
public:
    bool overlapped;
    int offset;
    int overlap_len;
    int diff;

    OverlapResult() : overlapped(false), offset(0), overlap_len(0), diff(0) {}
};

// 序列处理工具
class Sequence {
public:
    static string reverseComplement(string* seq) {
        string rc = *seq;
        reverse(rc.begin(), rc.end());
        for(int i = 0; i < rc.length(); i++) {
            switch(rc[i]) {
                case 'A': rc[i] = 'T'; break;
                case 'T': rc[i] = 'A'; break;
                case 'C': rc[i] = 'G'; break;
                case 'G': rc[i] = 'C'; break;
                default: break;
            }
        }
        return rc;
    }

    static char complement(char base) {
        switch(base) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
            case 'G': return 'C';
            default: return 'N';
        }
    }
};

// Overlap分析类
class OverlapAnalysis {
public:
    static OverlapResult analyze(Read* r1, Read* r2, int diffLimit = 5, int overlapRequire = 15, double diffPercentLimit = 0.2) {
        return analyze(r1->mSeq, r2->mSeq, diffLimit, overlapRequire, diffPercentLimit);
    }

    static OverlapResult analyze(string* r1, string* r2, int diffLimit, int overlapRequire, double diffPercentLimit) {
        string rcr2 = Sequence::reverseComplement(r2);
        int len1 = r1->length();
        int len2 = rcr2.length();
        const char* str1 = r1->c_str();
        const char* str2 = rcr2.c_str();

        int complete_compare_require = 50;
        int overlap_len = 0;
        int offset = 0;
        int diff = 0;

        // 正向搜索
        while (offset < len1 - overlapRequire) {
            overlap_len = min(len1 - offset, len2);
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

            diff = 0;
            int i = 0;
            for (i = 0; i < overlap_len; i++) {
                if (str1[offset + i] != str2[i]) {
                    diff += 1;
                    if (diff > overlapDiffLimit && i < complete_compare_require)
                        break;
                }
            }

            if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > complete_compare_require)) {
                OverlapResult ov;
                ov.overlapped = true;
                ov.offset = offset;
                ov.overlap_len = overlap_len;
                ov.diff = diff;
                return ov;
            }
            offset += 1;
        }

        // 反向搜索
        offset = 0;
        while (offset > -(len2 - overlapRequire)) {
            overlap_len = min(len1, len2 - abs(offset));
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

            diff = 0;
            int i = 0;
            for (i = 0; i < overlap_len; i++) {
                if (str1[i] != str2[-offset + i]) {
                    diff += 1;
                    if (diff > overlapDiffLimit && i < complete_compare_require)
                        break;
                }
            }

            if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > complete_compare_require)) {
                OverlapResult ov;
                ov.overlapped = true;
                ov.offset = offset;
                ov.overlap_len = overlap_len;
                ov.diff = diff;
                return ov;
            }
            offset -= 1;
        }

        OverlapResult ov;
        ov.overlapped = false;
        return ov;
    }
};

// 修改后的碱基校正类
class BaseCorrector {
public:
    struct CorrectionStats {
        long totalReads = 0;
        long totalPairs = 0;
        long totalBases = 0;
        long correctedReads = 0;
        long correctedBases = 0;

        void print() {
            cout << "=== 校正统计结果 ===" << endl;
            cout << "总reads数: " << totalReads << endl;
            cout << "总read pairs数: " << totalPairs << endl;
            cout << "总碱基数: " << totalBases << endl;
            cout << "校正的reads数: " << correctedReads << endl;
            cout << "校正的碱基数: " << correctedBases << endl;
            cout << "校正reads比例: " << (totalReads > 0 ? (double)correctedReads/totalReads*100 : 0) << "%" << endl;
            cout << "校正碱基比例: " << (totalBases > 0 ? (double)correctedBases/totalBases*100 : 0) << "%" << endl;
        }
    };

    static int correctByOverlapAnalysis(Read* r1, Read* r2, CorrectionStats& stats, int diffLimit = 5, int overlapRequire = 30, double diffPercentLimit = 0.2) {
        OverlapResult ov = OverlapAnalysis::analyze(r1, r2, diffLimit, overlapRequire, diffPercentLimit);
        return correctByOverlapAnalysis(r1, r2, stats, ov);
    }

    static int correctByOverlapAnalysis(Read* r1, Read* r2, CorrectionStats& stats, OverlapResult ov) {
        // 更新统计信息
        stats.totalReads += 2;  // R1 + R2
        stats.totalPairs += 1;
        stats.totalBases += r1->length() + r2->length();

        // 只处理有overlap的情况
        if (ov.diff == 0 || !ov.overlapped)
            return 0;

        int ol = ov.overlap_len;
        int start1 = max(0, ov.offset);
        int start2 = r2->length() - max(0, -ov.offset) - 1;

        const char* seq1 = r1->mSeq->c_str();
        const char* seq2 = r2->mSeq->c_str();

        int corrected = 0;
        bool r1Corrected = false;
        bool r2Corrected = false;

        for(int i = 0; i < ol; i++) {
            int p1 = start1 + i;
            int p2 = start2 - i;

            // 修改：忽略质量，只要碱基不一致就校正为N
            if(seq1[p1] != Sequence::complement(seq2[p2])) {
                // 将R1和R2对应位置都校正为N
                (*r1->mSeq)[p1] = 'N';
                (*r2->mSeq)[p2] = 'N';
                // 质量值设为最低
                (*r1->mQuality)[p1] = '!';  // Q0
                (*r2->mQuality)[p2] = '!';  // Q0

                corrected += 2;  // 校正了两个碱基
                r1Corrected = true;
                r2Corrected = true;
            }
        }

        // 更新校正统计
        if(corrected > 0) {
            if(r1Corrected && r2Corrected)
                stats.correctedReads += 2;
            else if(r1Corrected || r2Corrected)
                stats.correctedReads += 1;
            stats.correctedBases += corrected;
        }

        return corrected;
    }
};

// 简单的FASTQ读取器
class SimpleFastqReader {
private:
    ifstream file;

public:
    SimpleFastqReader(const string& filename) {
        file.open(filename);
        if(!file.is_open()) {
            cerr << "无法打开文件: " << filename << endl;
        }
    }

    ~SimpleFastqReader() {
        if(file.is_open()) file.close();
    }

    Read* read() {
        if(!file.is_open() || file.eof()) return nullptr;

        string line1, line2, line3, line4;
        if(!getline(file, line1) || !getline(file, line2) ||
           !getline(file, line3) || !getline(file, line4)) {
            return nullptr;
        }

        return new Read(new string(line1), new string(line2), new string(line4));
    }

    bool eof() { return file.eof(); }
};

// 主程序
int main(int argc, char* argv[]) {
    if(argc != 5) {
        cout << "用法: " << argv[0] << " <R1.fastq> <R2.fastq> <out_R1.fastq> <out_R2.fastq>" << endl;
        return 1;
    }

    string r1File = argv[1];
    string r2File = argv[2];
    string outR1File = argv[3];
    string outR2File = argv[4];

    SimpleFastqReader reader1(r1File);
    SimpleFastqReader reader2(r2File);
    ofstream writer1(outR1File);
    ofstream writer2(outR2File);

    if (!writer1.is_open() || !writer2.is_open()) {
        cerr << "错误：无法创建输出文件" << endl;
        return 1;
    }

    BaseCorrector::CorrectionStats stats;
    long processedPairs = 0;

    cout << "开始处理FASTQ文件..." << endl;
    cout << "输入文件: " << r1File << " 和 " << r2File << endl;
    cout << "输出文件: " << outR1File << " 和 " << outR2File << endl;

    while(!reader1.eof() && !reader2.eof()) {
        Read* r1 = reader1.read();
        Read* r2 = reader2.read();

        if(!r1 || !r2) break;

        // 进行overlap分析和校正
        BaseCorrector::correctByOverlapAnalysis(r1, r2, stats);

        // 写入校正后的结果
        writer1 << *r1->mName << "\n" << *r1->mSeq << "\n+\n" << *r1->mQuality << "\n";
        writer2 << *r2->mName << "\n" << *r2->mSeq << "\n+\n" << *r2->mQuality << "\n";

        // 清理内存
        delete r1->mName; delete r1->mSeq; delete r1->mQuality; delete r1;
        delete r2->mName; delete r2->mSeq; delete r2->mQuality; delete r2;

        processedPairs++;
        if (processedPairs % 100000 == 0) {
            cout << "已处理 " << processedPairs << " 对reads..." << endl;
        }
    }

    writer1.close();
    writer2.close();

    // 输出统计结果
    cout << "\n处理完成！" << endl;
    cout << "总处理read pairs: " << processedPairs << endl;
    stats.print();

    return 0;
}
