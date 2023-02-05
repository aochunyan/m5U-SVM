# m5U-SVM
# 1 Description
RNA 5-methyluridine (m5U) modifications are obtained by methylation at the C5 position of uridine catalyzed by pyrimidine methylation transferase, which is related to the development of human diseases. Accurate identification of m5U modification sites from RNA sequences can contribute to the understanding of their biological functions and the pathogenesis of related diseases. Despite the good performance of these computational methods, there are some drawbacks and limitations.
We have developed a novel predictor, m5U-SVM, based on multi-view features and machine learning algorithms to construct predictive models for identifying m5U modification sites from RNA sequences. Comparison results showed that the proposed model achieved more accurate and balanced performance than existing tools.

# 2. Availability
2.1. Webserver is available at: http://lab.malab.cn/~acy/m5USVM 

2.2 Datasets and source code are available at:
http://lab.malab.cn/~acy/m5USVM/Dataset.zip
https://github.com/aochunyan/m5U-SVM.git

# 3. Requirements
Before running, please make sure the following packages are installed in Python environment:
joblib==1.1.0 numpy==1.23.3 pandas==1.4.4 scikit-learn==1.1.2 gensim==3.4.0

# 4. Running
Changing working dir to m5USVM, and then running the following command:
python m5USVM.py -i test.fasta -m Full_transcript -o prediction_results.csv
-i: input file in fasta format
-m: sequence mode(Full_transcript or Mature_mRNA)
-o: output file name
