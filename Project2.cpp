#include <string>
#include <vector>
#include <algorithm>
#include <string>
#include "bio.h"
//Lisa Lipin 
//yay
using std::string; using std::vector; using std::reverse; using std::any_of;

bool IsValidDNASequence(const string & input) { //this checks if each character in the sequence is valid and returns true or false accordingly
    string possible_for_valid = "ACGT";
    for (char x : input) {
        if (possible_for_valid.find(x) == std::string::npos){
            return false;
        }
    }
    return true;
}

void GetReverseComplementSequence(const string & input,  string * const output) { //this voided function takes in a const reference to a string and 
// adds the opposite base to the output, a const pointer to a string, then reverses the output
    for (char x : input) {
        if (x == 'A') {
            *output += 'T';
        } else if (x == 'T') {
            *output += 'A';
        } else if (x == 'G') {
            *output += 'C';
        } else if (x == 'C') {
            *output += 'G';
        } 
    }
    reverse((*output).begin(), (*output).end()); //https://www.geeksforgeeks.org/reverse-a-string-in-c-cpp-different-methods/
}

string GetRNATranscript(const string & input) {// this changes every Thymine to a Uracil and returns a string
    string output = "";
    for (char x : input) {
        if (x == 'T') {
            output.push_back('U');
        } else {
            output.push_back(x);
        }
    }
    return output;
}

vector<vector<string>> GetReadingFramesAsCodons(const string & input) { //this takes a const reference to a string and returns a vector of six reading frames separated into codons. There are 3 offset each of the original and antiparallel sequence
    vector<vector<string>> reading_frames;
    vector <string> original_0, original_1, original_2, antiparallel_0, antiparallel_1, antiparallel_2;
    string codon = "";
    int offset = 0;
    string offset_string = "";
    string original_sequence = "";
    GetReverseComplementSequence(input, &original_sequence);
    string original = GetRNATranscript(original_sequence);
    while (offset < 3) { //for each offset, add on sets of 3. If the codon is too short, don't add it. This applies to antiparallel code as well
        offset_string = original.substr(offset);
        for (auto x : offset_string) {
            codon += x;
            if (static_cast<int>(codon.length()) == 3) {
                if (offset == 0) {
                    original_0.push_back(codon);
                } else if (offset == 1){
                    original_1.push_back(codon);
                } else if (offset == 2){
                    original_2.push_back(codon);
                }
                codon = "";
            }
        }
        offset++;
    }
    reading_frames.push_back(original_0);
    reading_frames.push_back(original_1);
    reading_frames.push_back(original_2);
    offset = 0;
    string reverse = GetRNATranscript(input);
    while (offset < 3) {
        offset_string = reverse.substr(offset);
        for (auto x:offset_string) {
            codon += x;
            if (static_cast<int>(codon.length()) == 3) {
                if (offset == 0) {
                    antiparallel_0.push_back(codon);
                } else if (offset == 1){
                    antiparallel_1.push_back(codon);
                } else if (offset == 2){
                    antiparallel_2.push_back(codon);
                }
                codon = "";
            }
        }
        offset++;
    }
    reading_frames.push_back(antiparallel_0);
    reading_frames.push_back(antiparallel_1);
    reading_frames.push_back(antiparallel_2);
    return reading_frames; 
}

string Translate(const vector<string> & codon_sequence) { //this converts each codon to the corresponding amino acid and returns a string which is the protein.
    string output = "";
    vector<string> A = {"GCU", "GCC", "GCA", "GCG"};
    vector<string> R = {"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"};
    vector<string> N = {"AUU", "AAC"};
    vector<string> D = {"GAU", "GAC"};
    vector<string> C = {"UGU", "UGC"};
    vector<string> Q = {"CAA", "CAG"};
    vector<string> E = {"GAA", "GAG"};
    vector<string> G = {"GGU", "GGC", "GGA", "GGG"};
    vector<string> H = {"CAU", "CAC"};
    vector<string> I = {"AUU", "AUC", "AUA"};
    vector<string> L = {"UUA", "UUG", "CUU","CUC", "CUA", "CUG"};
    vector<string> K = {"AAA", "AAG"};
    vector<string> M = {"AUG"};
    vector<string> F = {"UUU", "UUC"};
    vector<string> P = {"CCU", "CCC", "CCA", "CCG"};
    vector<string> S = {"UCU","UCC","UCA", "UCG", "AGU", "AGC"};
    vector<string> T = {"ACU", "ACC", "ACA", "ACG"};
    vector<string> W = {"UGG"};
    vector<string> Y = {"UAU", "UAC"};
    vector<string> V = {"GUU", "GUC", "GUA", "GUG"};
    vector<string> stop_codon = {"UAG", "UGA","UAA"};
    for (string x: codon_sequence) {
        if (any_of(A.begin(), A.end(),x)) { //https://www.tutorialspoint.com/cpp_standard_library/cpp_algorithm_any_of.htm
            output.push_back('A');
        } else if (any_of(R.begin(), R.end(),x)) {
            output.push_back('R');
        } else if (any_of(N.begin(), N.end(),x)) {
            output.push_back('N');
        } else if (any_of(D.begin(), D.end(),x)) {
            output.push_back('D');
        } else if (any_of(C.begin(), C.end(),x)) {
            output.push_back('C');
        } else if (any_of(Q.begin(), Q.end(),x)) {
            output.push_back('Q');
        } else if (any_of(E.begin(), E.end(),x)) {
            output.push_back('E');
        } else if (any_of(G.begin(), G.end(),x)) {
            output.push_back('G');
        } else if (any_of(H.begin(), H.end(),x)) {
            output.push_back('H');
        } else if (any_of(I.begin(), I.end(),x)) {
            output.push_back('I');
        } else if (any_of(L.begin(), L.end(),x)) {
            output.push_back('L');
        } else if (any_of(K.begin(), K.end(),x)) {
            output.push_back('K');
        } else if (any_of(M.begin(), M.end(),x)) {
            output.push_back('M');
        } else if (any_of(F.begin(), F.end(),x)) {
            output.push_back('F');
        } else if (any_of(P.begin(), P.end(),x)) {
            output.push_back('P');
        } else if (any_of(T.begin(), T.end(),x)) {
            output.push_back('S');
        } else if (any_of(T.begin(), T.end(),x)) {
            output.push_back('T');
        } else if (any_of(W.begin(), W.end(),x)) {
            output.push_back('W');
        } else if (any_of(Y.begin(), Y.end(),x)) {
            output.push_back('Y');
        } else if (any_of(V.begin(), V.end(),x)) {
            output.push_back('V');
        } else if (any_of(stop_codon.begin(), stop_codon.end(),x)) {
            output.push_back('*');
        }
    }
    return output;
}

string GetLongestOpenReadingFrame(const string & DNA_sequence) { //this checks each reading frame for valid proteins and returns the longest possible frame as a string
    string longest_reading_frame = "";
    vector<vector<string>> frames = GetReadingFramesAsCodons(DNA_sequence);
    for (vector<string> x:frames) {
        string frame = Translate(x);
        int i = 0;
        int old_i = 0;
        int stop = 0;
        bool continue_search = true;
        string open_frame = "";
        while(continue_search) {
            i = frame.find('M', old_i); 
            stop = frame.find('*', i);
            if (i == std::string::npos){ //if no starter codon is found or there is no stop codon following, move on to the next frame
                continue_search = false;
            } else if (stop == std::string::npos) {
                continue_search = false;
            } else {
                open_frame = frame.substr(old_i,i - old_i);
                if (open_frame.size() > longest_reading_frame.size()) {
                    longest_reading_frame = open_frame;
                }
            }
            old_i = i + 1;
            
    }
    return longest_reading_frame;
}
