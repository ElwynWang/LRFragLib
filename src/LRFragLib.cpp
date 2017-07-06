#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <time.h>
#include <valarray>
#include <vector>

#include "Parameters.h"

using namespace std;


typedef pair<string,double>PAIR;

int cmp(const PAIR& x , const PAIR& y)
{
    return x.second> y.second;
}

int cmp2(const PAIR& x , const PAIR& y)
{
    return x.second< y.second;
}

void ReadSeq(string &SeqName,string &AASeq, const string proname)
{
  string InputSequenceFile = string(InfileDic)+proname+".fasta.txt";
  ifstream infile(InputSequenceFile.c_str(),ios::in);
  string tmp("");
  while(getline(infile,tmp))
  {
    if(tmp[0] == '>')
    {
      SeqName = tmp.substr(1);
      while(1)
      {
        getline(infile,tmp);
        if(tmp[0] != '>' && !infile.eof())
          AASeq += tmp;
        else
          break;
      }
    }
    else if(infile.eof())
      break;
  }
  infile.close();
}

int PrintSeqInfo(const string &SeqName,const string &AASeq)
{
  cout << "******************** Protein sequence information ********************" << endl;
  cout << ">" << SeqName << endl << AASeq << endl;
  cout << "Sequence length = " << AASeq.length() << endl;
  return AASeq.length();
}

string getField(string str, string symbol, int Index)
{
    size_t pos=0;
    size_t len=str.length();
    size_t symbol_len=symbol.length();
    vector<int> findposes;

    if (symbol_len==0)
    {
        cout << "Wrong symbol while using getField function" << endl;
        exit(1);
    }

    while(pos<len)
    {
        int find_pos=str.find(symbol,pos);
        if (find_pos!=-1)
        {
            findposes.push_back(find_pos);
            pos=find_pos+symbol_len;
        }
        else 
            break;
    }

    size_t findnum=findposes.size();
    if (findnum!=0 and Index<=findnum)
    {
        if (Index==0)
            return str.substr(0,findposes[0]);
        else if (Index==findnum)
            return str.substr(findposes[Index-1]+1,len-findposes[Index-1]-1);
        else
            return str.substr(findposes[Index-1]+1,findposes[Index]-findposes[Index-1]-1);
    }
    else
    {
        cout << "Could not find symbol in the string" << endl;
        exit(1);
    }
}


class  RandomNumberGenerator{
public:
    RandomNumberGenerator();
    ~RandomNumberGenerator();
    void setup(int ij=1802, int kl=9373);
    double randomUniform();
    int randomInteger(int lower, int upper);
    double randomDouble(double lower, double upper);
private:
    vector<double> u_;
    double c_;
    double cd_;
    double cm_;
    int i97_;
    int j97_;
};

RandomNumberGenerator::RandomNumberGenerator()  
{
    setup();
}


RandomNumberGenerator::~RandomNumberGenerator()
{
}


void RandomNumberGenerator::setup(int ij, int kl)
{
    double s,t;
    int ii,i,j,k,l,jj,m;

    // allocate space for u_
    u_.resize(97);

    /*
    Handle the seed range errors
    First random number seed must be between 0 and 31328
    Second seed must have a value between 0 and 30081
    */
        
    if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081)
    {
        ij = 1802;
        kl = 9373;
    }

    i = (ij / 177) % 177 + 2;
    j = (ij % 177)       + 2;
    k = (kl / 169) % 178 + 1;
    l = (kl % 169);

    for (ii=0; ii<97; ii++) 
    {
        s = 0.0;
        t = 0.5;
        for (jj=0; jj<24; jj++) 
        {
            m = (((i * j) % 179) * k) % 179;
            i = j;
            j = k;
            k = m;
            l = (53 * l + 1) % 169;
            if (((l * m % 64)) >= 32)
                s += t;
            t *= 0.5;
        }
        u_[ii] = s;
    }

    c_    = 362436.0 / 16777216.0;
    cd_   = 7654321.0 / 16777216.0;
    cm_   = 16777213.0 / 16777216.0;
    i97_  = 97;
    j97_  = 33;
}

double RandomNumberGenerator::randomUniform()
{
    double uni;

    uni = u_[i97_-1] - u_[j97_-1];
    if (uni <= 0.0) uni++;

    u_[i97_-1] = uni;
    i97_--;

    if (i97_ == 0) i97_ = 97;
    j97_--;

    if (j97_ == 0) j97_ = 97;
    c_ -= cd_;

    if (c_ < 0.0) c_ += cm_;
    uni -= c_;

    if (uni < 0.0) uni++;

    return(uni);
}

int RandomNumberGenerator::randomInteger(int lower, int upper)
{
    return((int)(randomUniform() * (upper - lower + 1)) + lower);
}

double RandomNumberGenerator::randomDouble(double lower, double upper)  
{
    return((upper - lower) * randomUniform() + lower);
}




int main(int argc, char *argv[])
{ 
    if (argc!=2)
    {
        cout << "Usage: ./LRFragLib ProteinName\n";
        exit(1);
    }


    // Input Parameters
    const string proname=argv[1];

    int fragnum=0;
    
    string SeqName("");
    string AASeq("");

    ReadSeq(SeqName,AASeq,proname);
    const int samplelength=PrintSeqInfo(SeqName,AASeq);

    
    // Set Feature Calculation Parameters
    
    valarray<double>InitPhysicalArr(FactorNum);
    map<char,valarray<double> >InitPhysicalMap;
    
    for (int i=0;i<FactorNum;i++)
    {
        InitPhysicalArr[i]=InitFactorArr[i];
    }


    for(int i=0;i<AAIndexNum;i++)
    {
        slice InitPhysicalArrSlice (i*EachAAFactorNum,EachAAFactorNum,1);
        InitPhysicalMap.insert(pair<char,valarray<double> >(FactorArrAAIndex[i],InitPhysicalArr[InitPhysicalArrSlice]));
    }
       
    map< string,int >aaindex;
    for (int i=0;i<AAIndexNum;i++)
        aaindex.insert(pair<string,int>(string(FactorArrAAIndex).substr(i,1),i));


    vector<PAIR>pdbvecpair;
    vector<PAIR>pdblenvecpair;
    vector<PAIR>lenvecpair,tmplenvecpair;
    vector<PAIR>vecpair;
    vector<PAIR>rands;
    vector<PAIR>top100vecpair,top1000vecpair;
    set<string> duplicate; 
    set<string> HomolPdb;
    set<int> rmfrag; 

    string *lib7SSArr = new string [Lib7Num];
    string *lib8SSArr = new string [Lib8Num];
    string *lib9SSArr = new string [Lib9Num];
    string *lib10SSArr = new string [Lib10Num];
    string *DihedralAngs7Arr = new string [Lib7Num];
    string *DihedralAngs8Arr = new string [Lib8Num];
    string *DihedralAngs9Arr = new string [Lib9Num];
    string *DihedralAngs10Arr = new string [Lib10Num];


    RandomNumberGenerator *RNG = new RandomNumberGenerator();
    

    int *Top1pair = new int [samplelength-MinFragLen+1];  //only record libfragid
    string **Top15pair = new string* [samplelength-MinFragLen+1];
    for (int i=0;i<samplelength-MinFragLen+1;i++)
    {
        Top15pair[i] = new string [topcandidates_supply];
    }
                
    int EachposiNum[samplelength-MinFragLen+1];
    int eachlenposinum[samplelength-MinFragLen+1];
    for (int i=0;i<samplelength-MinFragLen+1;i++)
    {
        EachposiNum[i]=0;
        eachlenposinum[i]=0;
    }
    cout << "Set Parameters Done!" << endl; 
    if (rmhomol==true)
        cout << "Exclude Homologous: True!" << endl;
    else 
        cout << "Exclude Homologous: False!" << endl;

    
    // ***************************************************************  Begin  ************************************************** 
    
    for (int fraglength=MaxFragLen;fraglength>=MinFragLen;fraglength--)
    {   
        cout << "*************************  Calculation Begin *************************" << endl;
        cout << "Fragment length=" << fraglength << endl;
        // set common parameters
        double paras1[2*fraglength+2];
        double paras2[2*fraglength+1];
        double cutoff;
        
        char str_fraglen[10];
        sprintf(str_fraglen,"%d",fraglength);  

        // set parameters by length
        if (fraglength==7)
        {
            for (int newcnt=0;newcnt<2*fraglength+1;newcnt++)
            {
                paras1[newcnt]=paras7[newcnt];
                paras2[newcnt]=twoparas7[newcnt];
            }
            paras1[2*fraglength+1]=paras7[2*fraglength+1];
            cutoff=cutoff7;
            fragnum=Lib7Num;
        }

        else if (fraglength==8)
        {
            for (int newcnt=0;newcnt<2*fraglength+1;newcnt++)
            {
                paras1[newcnt]=paras8[newcnt];
                paras2[newcnt]=twoparas8[newcnt];
            }
            paras1[2*fraglength+1]=paras8[2*fraglength+1];
            cutoff=cutoff8;
            fragnum=Lib8Num;
        }

        else if (fraglength==9)
        {
            for (int newcnt=0;newcnt<2*fraglength+1;newcnt++)
            {
                paras1[newcnt]=paras9[newcnt];
                paras2[newcnt]=twoparas9[newcnt];
            }
            paras1[2*fraglength+1]=paras9[2*fraglength+1];
            cutoff=cutoff9;
            fragnum=Lib9Num;
        }

        else if (fraglength==10)
        {   
            for (int newcnt=0;newcnt<2*fraglength+1;newcnt++)
            {
                paras1[newcnt]=paras10[newcnt];
                paras2[newcnt]=twoparas10[newcnt];
            }
            paras1[2*fraglength+1]=paras10[2*fraglength+1];
            cutoff=cutoff10;
            fragnum=Lib10Num;
        }

        else
        {
            cout << "Error! Wrong Fragment Length\n";
            exit(1);
        }


        


        // **********************************************    Set Library Infomation  *****************************************************
        unsigned int fragcount=0;

        // Frag and Info
        string linefrag;
        string *FragmentArr = new string [fragnum];
        string *Frag2PdbArr = new string [fragnum];
        string *FragBeginArr = new string [fragnum];
        string dicfrag="./database/"+string(str_fraglen)+"-fraglist.txt";
        ifstream infilefrag(dicfrag.c_str(),ios::in);
        if (infilefrag.good())
        {
            while(getline(infilefrag,linefrag))
            {   
                *(FragmentArr+fragcount)= linefrag.substr(0,fraglength);
                *(Frag2PdbArr+fragcount)= linefrag.substr(fraglength+1,5);
                *(FragBeginArr+fragcount)= linefrag.substr(fraglength+7,linefrag.size()-fraglength-7);
                fragcount++;
            }
            infilefrag.close();
            fragcount=0;
        }
        else
        { 
            cout << "Could not open fragment library file!" << endl;
            exit(1);
        }


        // Slide frag
        string lineslide;
        string *SlidFragmentArr = new string [fragnum];
        string dicslide="./database/"+string(str_fraglen)+"-fragslidewindow7.txt";
        ifstream infileslide(dicslide.c_str(),ios::in);
        if (infileslide.good())
        {
            while(getline(infileslide,lineslide))
            {   
                *(SlidFragmentArr+fragcount)= lineslide;
                fragcount++;
            }
            infileslide.close();
            fragcount=0;
        }
        else
        {
            cout << "Could not open fragment slidewindow library!" << endl;
            exit(1);
        }



        // SS
        string liness;
        string *libSSArr= new string [fragnum];
        string dicss="./database/"+string(str_fraglen)+"-sslist.txt";
        ifstream infiless(dicss.c_str(),ios::in);
        if (infiless.good())
        {
            while(getline(infiless,liness))
            {           
                *(libSSArr+fragcount)= liness;

                if (fraglength==7)
                    *(lib7SSArr+fragcount)=liness;
                else if (fraglength==8)
                     *(lib8SSArr+fragcount)=liness;
                else if (fraglength==9)
                     *(lib9SSArr+fragcount)=liness;
                else 
                     *(lib10SSArr+fragcount)=liness;
                fragcount++;
            } 
            infiless.close();
            fragcount=0;
        }
        else
        {
            cout << "Could not open secondary structure library!" << endl;
            exit(1);
        }


        // Angle
        string lineang;
        double **LibAngleArr= new double* [fragnum];
        for (unsigned int i=0;i<fragnum;i++)
            LibAngleArr[i]= new double [fraglength*2];
        rmfrag.clear();
        string dicang="./database/"+string(str_fraglen)+"-anglelist.txt";
        ifstream infileang(dicang.c_str(),ios::in);
        if (infileang.good())
        {
            while(getline(infileang,lineang))
            {
                for (int j=0;j<fraglength;j++)
                {
                    LibAngleArr[fragcount][j*2]=atof(getField(lineang,",",3*j).c_str());
                    LibAngleArr[fragcount][j*2+1]=atof(getField(lineang,",",3*j+1).c_str());
                   
                    if (abs(atof(getField(lineang,",",3*j+2).c_str()))<160.00)
                        rmfrag.insert(fragcount);

                    if (fraglength==7)
                        *(DihedralAngs7Arr+fragcount)=lineang;
                    else if (fraglength==8)
                        *(DihedralAngs8Arr+fragcount)=lineang;
                    else if (fraglength==9)
                        *(DihedralAngs9Arr+fragcount)=lineang;
                    else 
                        *(DihedralAngs10Arr+fragcount)=lineang;
                }
                fragcount++;
            }
            infileang.close();
            fragcount=0;
        }
        else
        {
            cout << "Could not open rama angle library!" << endl;
            exit(1);
        }


        // **********************************************    Set sample infomation  *****************************************************
        // Set Parameters
        int pdbfragnum=samplelength-fraglength+1;
        int pairnum=pdbfragnum*fragnum; 
        const int LongFragmentLength=fraglength;


        // SPFrag and SPSlideFrag
        string linesppep;
        string *SampleFragmentArr = new string [pdbfragnum];
        string *SampleSlidFragmentArr = new string [pdbfragnum];
        string dicsppep=string(InfileDic)+proname+".fasta.txt";  
        ifstream infilesppep(dicsppep.c_str(),ios::in);
        string Sequence("");
        if (infilesppep.good())
        {
            while(getline(infilesppep,linesppep))
            {   
                if (linesppep[0]!='>')
                    Sequence+=linesppep;

            }
            infilesppep.close();
        }
        else
        {
            cout << "Could not open sample fasta file!" << endl;
            exit(1);
        }
   
        for (int seqindex=0;seqindex<samplelength-fraglength+1;seqindex++)
        {
            string frag_(fraglength,'A');
           
            string frag=frag_;
            for (int j=0;j<fraglength;j++)
                frag[j]=Sequence[seqindex+j];
            *(SampleFragmentArr+fragcount)= frag;
            fragcount++;
        }    
        fragcount=0;


        string SlideSequence=Sequence.substr(0,1)+Sequence.substr(0,1)+Sequence.substr(0,1)+Sequence+Sequence.substr(samplelength-1,1)+Sequence.substr(samplelength-1,1)+Sequence.substr(samplelength-1,1);

        for (int seqindex=0;seqindex<samplelength-fraglength+1;seqindex++)
        {
            string frag_(fraglength+6,'A');
            string frag=frag_;
            for (int j=0;j<fraglength+6;j++)
                frag[j]=SlideSequence[seqindex+j];          
            *(SampleSlidFragmentArr+fragcount)= frag;
            fragcount++;
        }    
        fragcount=0;


        // SPSS
        string linespss; 
        string *SampleSSArr= new string [pdbfragnum];
        string dicspss=string(InfileDic)+proname+".fasta.ss";  
        ifstream infilespss(dicspss.c_str(),ios::in);
        int sscount=0; 
        string sst(fraglength,'C');
        string SST = sst;
        string SecondSequence="";
        if (infilespss.good())
        {
            while(getline(infilespss,linespss))
            {
                if (linespss[7]=='H' or linespss[7]=='G' or linespss[7]=='I'or linespss[7]=='B'or linespss[7]=='E' or linespss[7]=='C'or linespss[7]=='S' or linespss[7]=='T')
                    SecondSequence+=linespss.substr(7,1);
                else
                {
                    cout << " Could not recognize secondary structure symbol!" << endl;
                    exit(1);
                }

            }

            infilespss.close();
        }
        else
        {
            cout << "Could not open sample secondary structure file!" << endl;
            exit(1);
        }

        for (int seqindex=0;seqindex<samplelength-fraglength+1;seqindex++)
        {
            for(int j=0;j<fraglength;j++)
            {
                SST[j]=SecondSequence[seqindex+j];
                if(SST[j]=='G' or SST[j]=='I')
                    SST[j]='H';
                if(SST[j]=='B')
                    SST[j]='E';
                if(SST[j]=='S')
                    SST[j]='C';
            }
            *(SampleSSArr+sscount)=SST; 
            sscount++;
        }


        // SpAngle
        string linespang;
        double *SampleAngleArr= new double [samplelength*2];
        string dicspang=string(InfileDic)+proname+".fasta.phipsi";
        ifstream infilespang(dicspang.c_str(),ios::in);
        int spanglinecnt=0;
        while(getline(infilespang,linespang))
        {
            if (spanglinecnt>0)
            { 
                SampleAngleArr[(spanglinecnt-1)*2]=atof(getField(linespang, "\t", 4).c_str());
                SampleAngleArr[(spanglinecnt-1)*2+1]=atof(getField(linespang, "\t", 5).c_str());
                //SampleAngleArr[(spanglinecnt-1)*2]=atof(linespang.substr(17,6).c_str());  //old version of spine-x
                //SampleAngleArr[(spanglinecnt-1)*2+1]=atof(linespang.substr(25,6).c_str());  //old version of spine-x
            }
            spanglinecnt++;
        }
        infilespang.close();


        // Remove Homologous
        if (rmhomol)       
        {
            string dichomol=string(InfileDic)+proname+".homol";
            ifstream infilehomol(dichomol.c_str(),ios::in);
            string linehomol;
            int linecnt=0;
            HomolPdb.clear();
            
            while(getline(infilehomol,linehomol))
            {
                linecnt++;
                if (linecnt==25 and linehomol.substr(1,26)=="***** No hits found ******")
                    break;
                if (linecnt>=28 and linehomol=="")
                    break;
                if (linecnt>=28)
                    HomolPdb.insert(linehomol.substr(0,5));
            }
        }



        /*****************************************************       Begin to do calculation     ***********************************/
        for (int pdbfragcnt=0;pdbfragcnt<pdbfragnum;pdbfragcnt++)
        {
            string longfragment1=*(SampleSlidFragmentArr+pdbfragcnt); 

            string spss=*(SampleSSArr+pdbfragcnt);
            bool docalculation=true;
            int hcnt,ecnt,ccnt;
            hcnt=ecnt=ccnt=0;
            for (int ssidx=0;ssidx<fraglength;ssidx++)
            {
                if (spss[ssidx]=='H')
                    hcnt+=1;
                else if (spss[ssidx]=='E')
                    ecnt+=1;
                else if (spss[ssidx]=='C' or spss[ssidx]=='T')
                    ccnt+=1;
                else 
                {
                    cout << "Could not recognize secondary structure symbol!" << endl;
                    exit(1);
                }
            }

            
            if (fraglength>=9 and ccnt<5)             // helix:9,8,7 strand:7 coil:10,9,8,7
                docalculation=false;
            else if (fraglength==8 and ecnt>=5)
                docalculation=false;

            if (!docalculation)
                continue;
            /*
            if (ccnt>=cminnum)     omit other for len=8,7             
                docalculation=true;
            else if (hcnt>=hminnum and fraglength<=hmaxlen)
                docalculation=true;
            else if (ecnt>=eminnum and fraglength<=emaxlen)
                docalculation=true;

            if (!docalculation)
                continue;
            */

            #pragma omp parallel for
            for (int libfragcnt=0;libfragcnt<fragnum;libfragcnt++)
            {
                bool skipfrag=false;                                              
                if (rmhomol)
                {
                    string frag2pdb=Frag2PdbArr[libfragcnt];
                    set<string>::iterator setiter;
                    if ( HomolPdb.size() == 0 )
                        skipfrag=false; 
                    else if ( (setiter = HomolPdb.find(frag2pdb)) != HomolPdb.end())
                    {
                        skipfrag=true;
                    }
                }
                set<int>::iterator setiterfrag;
                if ( (setiterfrag = rmfrag.find(libfragcnt)) != rmfrag.end())
                {
                    skipfrag=true;
                }

                if (skipfrag)
                    continue;                                              

                string longfragment2=*(SlidFragmentArr+libfragcnt);
                double features[2*fraglength+1];

                if (longfragment1.size()==fraglength+6 and longfragment2.size()==fraglength+6)
                {
        
                    //**********************************************Part 1: Physical Factor*************************************************
                    if (ppflag)
                    {
                        for(int j=(ShortFragmentLength-1)/2;j<fraglength+(ShortFragmentLength-1)/2;j++)
                        {
                            string fragment1=longfragment1.substr(j-(ShortFragmentLength-1)/2,j+(ShortFragmentLength-1)/2+1);
                            string fragment2=longfragment2.substr(j-(ShortFragmentLength-1)/2,j+(ShortFragmentLength-1)/2+1);
           
                            valarray<double>Arr1(10*ShortFragmentLength),Arr2(10*ShortFragmentLength),tmp1(10),tmp2(10);   // 1 is for target fragment and 2 is for the library fragment
                     
                            map< char,valarray<double> >::iterator map_it;
                            double sum1,sum2;
                            double numerator,denominator1,denominator2,r;

                            // build arr1 and arr2 according InitPhysicalMap
                            for(int i=0;i<ShortFragmentLength;i++)
                            {
                                for (map_it=InitPhysicalMap.begin();map_it!=InitPhysicalMap.end();map_it++)
                                {
                                    if(fragment1[i]==map_it->first)
                                    {
                                        tmp1=map_it->second;
                                        for(int k=0;k<10;k++)
                                            Arr1[i*10+k]=tmp1[k];
                                    }
                                    if(fragment2[i]==map_it->first)
                                    {
                                        tmp2=map_it->second;
                                        for(int k=0;k<10;k++)
                                            Arr2[i*10+k]=tmp2[k];
                                    }
                                }
                            }
                       
                            sum1=Arr1.sum();
                            sum2=Arr2.sum();
               
                            Arr1=Arr1-sum1/(10*ShortFragmentLength);
                            Arr2=Arr2-sum2/(10*ShortFragmentLength);

                            numerator=(Arr1*Arr2).sum();
                            denominator1=(Arr1*Arr1).sum();
                            denominator2=(Arr2*Arr2).sum();
                       
                            r=(numerator/(sqrt(denominator1*denominator2))+1.00)/2.0;
                            features[j-(ShortFragmentLength-1)/2]=r;
                        }
                    }
                    else
                    {
                        for(int j=(ShortFragmentLength-1)/2;j<fraglength+(ShortFragmentLength-1)/2;j++)
                            features[j-(ShortFragmentLength-1)/2]=0.0;
                    }

                    //***************************************Part 2: Primary Sequence***************************************************
                    if(psflag)
                    {   
                        for(int j=(ShortFragmentLength-1)/2;j<fraglength+(ShortFragmentLength-1)/2;j++)
                        {
                            string fragment1=longfragment1.substr(j-(ShortFragmentLength-1)/2,j+(ShortFragmentLength-1)/2+1);
                            string fragment2=longfragment2.substr(j-(ShortFragmentLength-1)/2,j+(ShortFragmentLength-1)/2+1);
                            int sum=0;
            
                            for(int k=0;k<ShortFragmentLength;k++)
                            {
                                map<string,int>::iterator it1,it2;
                                it1=aaindex.find(fragment1.substr(k,1));
                                it2=aaindex.find(fragment2.substr(k,1));
                                int p=it1->second;
                                int q=it2->second;
                                sum+=Blosum62[p][q];                    
                            }
                    
                            double normsum=(double(sum)+double(4*ShortFragmentLength+1))/double(15*ShortFragmentLength+1);
                            int posi=fraglength+j-(ShortFragmentLength-1)/2;
                            features[posi]=normsum;
                        }
                    }
                    else
                    {
                        for(int j=(ShortFragmentLength-1)/2;j<fraglength+(ShortFragmentLength-1)/2;j++)
                            features[fraglength+j-(ShortFragmentLength-1)/2]=0.0;
                    }
                }
                else
                {
                    cout << "Error: Program exited!!! Wrong sequence length for calculation!!!"  << endl;
                    exit(1);
                }
                //***************************************Part 3: Secondary Structure ***************************************************
                if (ssflag)
                {
                    string SS1=*(SampleSSArr+pdbfragcnt);
                    string SS2=*(libSSArr+libfragcnt);

                    int SimilarityScore = 0.0;
                    for(int k=0;k<fraglength;k++)
                    {
                        if(SS1[k] == SS2[k])
                            SimilarityScore += 3.0;
                        else if(SS1[k] == 'C' or SS2[k] == 'C' or SS2[k] == 'T' or SS2[k] == 'T')
                            SimilarityScore += 1.0;
                        else if(SS1[k] == 'H' or SS1[k] == 'E' or SS2[k] == 'H' or SS2[k] == 'E')
                            SimilarityScore += 0.0;
                        else
                        {    
                            cerr << "Error: Program exited!!! Can't find the corresponding secondary structure element!!!" << endl;
                            exit(1);
                        }
                    }

                    double normolscore=SimilarityScore/(3.000*fraglength);
                    features[2*fraglength]=normolscore;
                }
                else
                {
                    features[2*fraglength]=0.0;
                }


                double value1=paras1[0];
                double value2=paras2[0];
                for (int cnt=0;cnt<2*fraglength;cnt++)
                {   
                    if (cnt<fraglength)
                    {
                        value1+=features[cnt]*paras1[cnt+1];
                        value2+=features[cnt]*paras2[cnt+1];
                    }
                    else
                    {
                        value1+=features[cnt]*paras1[cnt+1];
                        value2+=features[cnt]*paras2[cnt+1];
                    }            
                }


                value1+=features[2*fraglength]*paras1[2*fraglength+1];
                double trueval1=1.000/(exp(-value1)+1.000);
                double trueval2=1.000/(exp(-value2)+1.000);
                
                char str_val2[100];
                sprintf(str_val2,"%f",trueval2);
                char str_pdbfragcnt[100];
                sprintf(str_pdbfragcnt, "%d", pdbfragcnt);
                char str_libfragcnt[100];
                sprintf(str_libfragcnt, "%d", libfragcnt);

                if (trueval1>cutoff)
                {
                    string firstitem=string(str_pdbfragcnt)+","+string(str_libfragcnt)+","+string(str_val2);
                    #pragma omp critical
                    vecpair.push_back(make_pair(firstitem,trueval1));
                }
            }

            //Sort all fragments of each position for each length
            sort(vecpair.begin(),vecpair.end(),cmp); 

            // Prepare for supplement step
            if (vecpair.size()==0)
            {
                cout << "Could not find any fragment!!!" << endl;
                exit(1);
            }

            if (fraglength==7 and vecpair.size()!=0)       
            {   
                if (vecpair.size()<candidates_supply)
                {
                    for (int i=0;i<vecpair.size();i++)
                    {
                        int spidx=atoi(getField(vecpair[i].first, ",", 0).c_str());
                        int libidx=atoi(getField(vecpair[i].first, ",", 1).c_str());

                        double anglescore=0.0;
                        for (int m=0;m<fraglength;m++)
                        {
                            double deltaphi=min(abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]),360.0-abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]));
                            double deltapsi=min(abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]),360.0-abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]));
                            anglescore=anglescore+deltaphi+deltapsi;
                        }
                        char str_vecsecond[100];
                        sprintf(str_vecsecond,"%f",vecpair[i].second);
                        top100vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                        top1000vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                    }
                }
                else if (vecpair.size()<candidates_enrich)
                {
                    for (int i=0;i<vecpair.size();i++)
                    {
                        int spidx=atoi(getField(vecpair[i].first, ",", 0).c_str());
                        int libidx=atoi(getField(vecpair[i].first, ",", 1).c_str());
                        
                        double anglescore=0.0;
                        
                        char str_vecsecond[100];
                        sprintf(str_vecsecond,"%f",vecpair[i].second);

                        for (int m=0;m<fraglength;m++)
                        {
                            double deltaphi=min(abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]),360.0-abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]));
                            double deltapsi=min(abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]),360.0-abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]));
                            anglescore=anglescore+deltaphi+deltapsi;
                        }
                        if (i<candidates_supply)
                        {
                            top100vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                            top1000vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                        }
                        else
                            top1000vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                        
                    }
                }
                else
                {
                    for (int i=0;i<candidates_enrich;i++)
                    {
                        int spidx=atoi(getField(vecpair[i].first, ",", 0).c_str());
                        int libidx=atoi(getField(vecpair[i].first, ",", 1).c_str());
                       
                        double anglescore=0.0;

                        char str_vecsecond[100];
                        sprintf(str_vecsecond,"%f",vecpair[i].second);

                        for (int m=0;m<fraglength;m++)
                        {
                            double deltaphi=min(abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]),360.0-abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]));
                            double deltapsi=min(abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]),360.0-abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]));
                            anglescore=anglescore+deltaphi+deltapsi;
                        }
                        if (i<candidates_supply)
                        {
                            top100vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                            top1000vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                        }
                        else
                            top1000vecpair.push_back(make_pair(string(vecpair[i].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","+FragBeginArr[libidx],anglescore));
                        
                    }
                }

                // get anglescore top15 from top100 for each position of each length
                int minnum;
                if (top100vecpair.size()<topcandidates_supply)
                    minnum=top100vecpair.size();
                else
                    minnum=topcandidates_supply;
                partial_sort(top100vecpair.begin(),top100vecpair.begin()+minnum,top100vecpair.end(),cmp2);  
                
                int spidx=atoi(getField(top100vecpair[0].first, ",", 0).c_str());
               
                
                for  (int i=0;i<minnum;i++)
                {
                    char str_vecsecond[100];
                    sprintf(str_vecsecond,"%f",top100vecpair[i].second);
                    Top15pair[spidx][i]=string(top100vecpair[i].first)+","+string(str_vecsecond);
                }
                  
                vector<PAIR>().swap(top100vecpair);

            }



            // Choose fragments step1 (candition1) --level of each position of a peptide
            #pragma omp parallel for 
            for (int l=0;l<vecpair.size();l++)   
            {
                int spidx=atoi(getField(vecpair[l].first, ",", 0).c_str());
                int libidx=atoi(getField(vecpair[l].first, ",", 1).c_str());
               
                double anglescore=0.0;
                double ratio=double(l+1.00)/(vecpair.size());
                
                char str_ratio[100];
                sprintf(str_ratio,"%f",ratio);
                char str_vecsecond[100];
                sprintf(str_vecsecond,"%f",vecpair[l].second);
                
                for (int m=0;m<fraglength;m++)
                {
                    double deltaphi=min(abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]),360.0-abs(SampleAngleArr[spidx*2+m*2]-LibAngleArr[libidx][m*2]));
                    double deltapsi=min(abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]),360.0-abs(SampleAngleArr[spidx*2+m*2+1]-LibAngleArr[libidx][m*2+1]));
                    anglescore=anglescore+deltaphi+deltapsi;
                }
      

                #pragma omp critical
                tmplenvecpair.push_back(make_pair(string(vecpair[l].first)+","+string(str_vecsecond)+","+FragmentArr[libidx]+","+Frag2PdbArr[libidx]+","\
                                                   +FragBeginArr[libidx]+","+string(str_ratio),anglescore));  // ratio need to be removed at the end
            }    
 
            if (vecpair.size()!=0)
            {
                vector<PAIR>().swap(vecpair); 
                lenvecpair.insert(lenvecpair.end(),tmplenvecpair.begin(),tmplenvecpair.end()); 
            }

            vector<PAIR>().swap(tmplenvecpair);
        }

        //Free of memory
        
        for (int i=0;i<fragnum;i++)
            delete [] LibAngleArr[i];

        delete[] SlidFragmentArr;
        delete[] libSSArr;
        delete[] LibAngleArr;
        delete[] FragmentArr;
        delete[] Frag2PdbArr;
        delete[] FragBeginArr;
        FragmentArr=NULL;
        Frag2PdbArr=NULL;
        FragBeginArr=NULL;
        SlidFragmentArr=NULL;
        libSSArr=NULL;
        LibAngleArr=NULL;


        // Choose fragments step2-- level of all positions in a peptide
        for (int cycle=0;cycle<pdbfragnum*2;cycle++)
        {  
            int lensize=lenvecpair.size();
            int *idxrand = new int [lensize];
            for (int i =0;i<lensize;i++)
                idxrand[i]=1;
            int randcnt=0;
            int selectnum=min(EachRoundNum,lensize/10);
            
            while(randcnt<selectnum)
            {
                int randomidx=RNG->randomInteger(0,lensize-1);

                if (idxrand[randomidx]==1)
                {
                    randcnt++;
                    idxrand[randomidx]=0;
                    if (atof(getField(lenvecpair[randomidx].first, ",", 7).c_str())<ThreeFactorAccRatio and atof(getField(lenvecpair[randomidx].first, ",", 2).c_str())>TwoFactoreAccVal)
                    {
                        size_t pos7=lenvecpair[randomidx].first.rfind(",");
                        rands.push_back(make_pair(lenvecpair[randomidx].first.substr(0,pos7), lenvecpair[randomidx].second));
                    }
                    else if (atof(getField(lenvecpair[randomidx].first, ",", 3).c_str())>ThreeFactorAccVal and atof(getField(lenvecpair[randomidx].first, ",", 2).c_str())>TwoFactoreAccVal)
                    {
                        size_t pos7=lenvecpair[randomidx].first.rfind(",");
                        rands.push_back(make_pair(lenvecpair[randomidx].first.substr(0,pos7), lenvecpair[randomidx].second));
                    }
                } 
            }

            // Select top 10
            if (rands.size()>=10)
            {
                partial_sort(rands.begin(),rands.begin()+10,rands.end(),cmp2);
                for (int i=0;i<10;i++)
                {
                    pdblenvecpair.push_back(make_pair(rands[i].first,rands[i].second));
                    int spidx=atoi(getField(rands[i].first, ",", 0).c_str());

                    eachlenposinum[spidx]++;
                }
            }
            else
            {
                for (int i=0;i<rands.size();i++)
                {
                    pdblenvecpair.push_back(make_pair(rands[i].first,rands[i].second));
                    int spidx=atoi(getField(rands[i].first, ",", 0).c_str());
                    
                    eachlenposinum[spidx]++;
                }
            }
            
            vector<PAIR>().swap(rands);
            delete [] idxrand;
            idxrand=NULL; 
        }
        vector<PAIR>().swap(lenvecpair);

        if (fraglength==7)
        {
            double Anglescorecompare[samplelength-6];
            for (int i=0;i<samplelength-6;i++)
                Anglescorecompare[i]=10000.00; 

            for (int i=0;i<pdblenvecpair.size();i++)
            {
                int spidx=atoi(getField(pdblenvecpair[i].first, ",", 0).c_str());
                int libidx=atoi(getField(pdblenvecpair[i].first, ",", 1).c_str());
                
                double angscore=pdblenvecpair[i].second;
                if (angscore<Anglescorecompare[spidx])
                {
                   Top1pair[spidx]=libidx;   // only record library index  
                   Anglescorecompare[spidx]=angscore; 
                }
            }
            for (int i=0;i<samplelength-6;i++)
            {
                if (Anglescorecompare[i]==10000.00)
                    Top1pair[i]=-1;
            }
        }


        // Integrate chosen fragments for each length
        for (int i=0;i<pdblenvecpair.size();i++)
        {
            int spidx=atoi(getField(pdblenvecpair[i].first, ",", 0).c_str());
            
            char str_ang[100];
            sprintf(str_ang,"%f",pdblenvecpair[i].second);
            
            // Length=10,9,8
            if (fraglength!=MinFragLen and EachposiNum[spidx]<=MaxPosSelcNum)
            {
                
                if (pdblenvecpair[i].second<AngleCriterion and eachlenposinum[spidx]>TakeallMinCandiNum)
                {
                    pdbvecpair.push_back(make_pair(string(pdblenvecpair[i].first)+","+string(str_ang)+","+string(str_fraglen),double(spidx)));
                    EachposiNum[spidx]++;
                }
                else if (pdblenvecpair[i].second<AngleCriterion and eachlenposinum[spidx]>RandSelecMinCandiNum)
                {
                    double random1=RNG->randomDouble(0,1);
                    double random2=RNG->randomDouble(0,1);
                    
                    if (random1>random2)
                    {
                        pdbvecpair.push_back(make_pair(string(pdblenvecpair[i].first)+","+string(str_ang)+","+string(str_fraglen),double(spidx)));
                        EachposiNum[spidx]++;
                    }
                }
            }
            // Length=7
            else if (fraglength==MinFragLen and EachposiNum[spidx]<=MaxPosSelcNum)
            {
                pdbvecpair.push_back(make_pair(string(pdblenvecpair[i].first)+","+string(str_ang)+","+string(str_fraglen),double(spidx)));
                EachposiNum[spidx]++;
                string dupidx=getField(pdblenvecpair[i].first, ",", 0)+","+getField(pdblenvecpair[i].first, ",", 1);
               
                duplicate.insert(dupidx);  //record spidx and libidx
            }           
        }
        vector<PAIR>().swap(pdblenvecpair);
        memset(eachlenposinum,samplelength-MinFragLen+1,0);

        delete[] SampleFragmentArr;
        delete[] SampleSlidFragmentArr;
        delete[] SampleSSArr;
        delete[] SampleAngleArr;

        SampleFragmentArr=NULL;
        SampleSlidFragmentArr=NULL;
        SampleSSArr=NULL;
        SampleAngleArr=NULL;
        

        cout << "Done! Length=" << fraglength << endl;
    }



    // Enrichment Step
    // Drms
    int fragcount=0;
    string linedrms; 
    string *DrmsArr = new string [DrmsItemNum];
    string dicdrms="./database/drmslist_minlen.txt";
    ifstream infiledrms(dicdrms.c_str(),ios::in);
    if (infiledrms.good())
    {
        while(getline(infiledrms,linedrms))
        {
            *(DrmsArr+fragcount)= linedrms;
            fragcount++;
        }
        infiledrms.close();
    }
    else
    {
        cout << "Could not open library drms file!" << endl; 
        exit(1);
    }
    
    #pragma omp parallel for 
    for (int i=0;i<top1000vecpair.size();i++)
    {
        string dupidx=getField(top1000vecpair[i].first, ",", 0)+","+getField(top1000vecpair[i].first, ",", 1);
      
        set<string>::iterator dupiter;
        dupiter = duplicate.find(dupidx);
        int spidx=atoi(getField(top1000vecpair[i].first, ",", 0).c_str());
     
        if ( dupiter == duplicate.end() and Top1pair[spidx]!=-1)
        {
            int frag1idx=Top1pair[spidx];
            int frag2idx=atoi(getField(top1000vecpair[i].first, ",", 1).c_str());
          
            valarray<double>Arr1(21),Arr2(21),Arr3(21),Arr4(21);
            
            char str_ang[100];
            sprintf(str_ang,"%f",top1000vecpair[i].second);

            for (int disnum=0;disnum<21;disnum++)
            {
                int pos1=frag1idx*21+disnum;
                int pos2=frag2idx*21+disnum;
                double dis1=atof((*(DrmsArr+pos1)).c_str());
                double dis2=atof((*(DrmsArr+pos2)).c_str());
                Arr1[disnum]=dis1;
                Arr2[disnum]=dis2;
            }

            Arr3=Arr1-Arr2;
            Arr4=Arr3*Arr3;
            double dist=sqrt(Arr4.sum()/21.0);

            if (dist<0.35)  
            {
                #pragma omp critical 
                {
                    pdbvecpair.push_back(make_pair(string(top1000vecpair[i].first)+","+string(str_ang)+",7",\
                                        atof(getField(top1000vecpair[i].first, ",", 0).c_str())));
                    EachposiNum[spidx]++;
                }               
            }
        }                             
    }


    vector<PAIR>().swap(top1000vecpair);
    delete [] DrmsArr;
    DrmsArr=NULL;


    // Supplement Step
    for (int i=0;i<samplelength-MinFragLen+1;i++)
    {
        int cnt=0;
        while (EachposiNum[i]<MinPosSelcNum and Top15pair[i][cnt]!="" and cnt<MinPosSelcNum)
        {
            string dupidx=getField(Top15pair[i][cnt], ",", 0)+","+getField(Top15pair[i][cnt], ",", 1);
           
            set<string>::iterator dupiter;
            dupiter=duplicate.find(dupidx);
            if ( dupiter == duplicate.end())
            {
                pdbvecpair.push_back(make_pair(Top15pair[i][cnt]+",7",atof(getField(Top15pair[i][cnt], ",", 0).c_str())));
               
                EachposiNum[i]++;
            }
            cnt++;   
        }
    }


    // Output result
    string finallibdic=string(OutfileDic)+proname+".lib";
    ofstream finallib(finallibdic.c_str());
    string finalfragdic=string(OutfileDic)+proname+".frag";
    ofstream finalfrag(finalfragdic.c_str());

    finallib << "#Query_pos Libindex Fragment PDBid Fragment_pos Score1 Score2 Anglescore FragmentLength" << endl;

    sort(pdbvecpair.begin(),pdbvecpair.end(),cmp2);

    for (int i=0;i<pdbvecpair.size();i++)
    {
        finallib << int(pdbvecpair[i].second) << "\t" << getField(pdbvecpair[i].first,",",1) << "\t" \
                 << getField(pdbvecpair[i].first,",",4) << "\t" << getField(pdbvecpair[i].first,",",5) << "\t" \
                 << getField(pdbvecpair[i].first,",",6) << "\t" << getField(pdbvecpair[i].first,",",2) << "\t" \
                 << getField(pdbvecpair[i].first,",",3) << "\t" << getField(pdbvecpair[i].first,",",7) << "\t" \
                 << getField(pdbvecpair[i].first,",",8) << endl;

    }

    int cntfrag=0;
    finalfrag << "position:          0" << " neighbors:          " << EachposiNum[0] << endl << endl;
    for (int i=0;i<pdbvecpair.size();i++)
    {
        int spidx=atoi(getField(pdbvecpair[i].first, ",", 0).c_str());
        if (spidx-1==cntfrag)
        {
            finalfrag << "position:          " << spidx << " neighbors:          " << EachposiNum[spidx] << endl << endl;
            cntfrag++;
        }
        
        int length=atoi(getField(pdbvecpair[i].first, ",", 8).c_str());
        int libidx=atoi(getField(pdbvecpair[i].first, ",", 1).c_str());
      
        string ss,angles;
        if (length==7)
        {
            ss=lib7SSArr[libidx];
            angles=DihedralAngs7Arr[libidx];
        }
        else if (length==8)
        {
            ss=lib8SSArr[libidx];
            angles=DihedralAngs8Arr[libidx];
        }
        else if (length==9)
        {
            ss=lib9SSArr[libidx];
            angles=DihedralAngs9Arr[libidx];
        }
        else 
        {
            ss=lib10SSArr[libidx];
            angles=DihedralAngs10Arr[libidx];
        }

        for (int k=0;k<length;k++)
        {
            finalfrag << " " << getField(pdbvecpair[i].first,",",5).substr(0,4) << " " << getField(pdbvecpair[i].first,",",5)[4] \
                        << "  " << setw(4) << atoi(getField(pdbvecpair[i].first,",",6).c_str())+k << " " << getField(pdbvecpair[i].first,",",4)[k]\
                        << " " << ss[k] << " " << setw(9) << setprecision(3) << setiosflags(ios::fixed) << atof(getField(angles,",",k*3+0).c_str()) \
                        << setw(9) << setprecision(3) << setiosflags(ios::fixed) << atof(getField(angles,",",k*3+1).c_str()) \
                        << setw(9) << setprecision(3) << setiosflags(ios::fixed) << atof(getField(angles,",",k*3+2).c_str()) << endl;  
        }
        finalfrag << endl;
    }



    // free of memory
    vector<PAIR>().swap(pdbvecpair);
    
    for (int i=0;i<samplelength-6;i++)
        delete [] Top15pair[i];
    
    delete [] Top15pair;
    delete [] Top1pair;
    Top15pair=NULL;
    Top1pair=NULL;
    
    delete [] lib7SSArr;
    delete [] lib8SSArr;
    delete [] lib9SSArr;
    delete [] lib10SSArr;
    delete [] DihedralAngs7Arr;
    delete [] DihedralAngs8Arr;
    delete [] DihedralAngs9Arr;
    delete [] DihedralAngs10Arr;
    lib7SSArr=NULL;
    lib8SSArr=NULL;
    lib9SSArr=NULL;
    lib10SSArr=NULL;
    DihedralAngs7Arr=NULL;
    DihedralAngs8Arr=NULL;
    DihedralAngs9Arr=NULL;
    DihedralAngs10Arr=NULL;
    
    cout << "*******************************  Done!  ******************************" << endl;
    
    return 0;
}    












