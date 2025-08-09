# 为DNNGP制作的输入格式转换工具 Input format conversion tool for DNNGP
DNNG是一个用于基因组预测的深度神经网络模型
原项目位置： > [DNNGP: Deep neural network for genomic prediction](https://github.com/AIBreeding/DNNGP)

本项目的作用是将hmp.txt文件和csv文件转化为输入模型需要的pkl以及tsv文件，并包含了对齐样本的功能。

请根据流程图内的使用顺序使用。各类程序内有其具体使用说明。

下面是各类文件或文件夹的使用说明
- ### SNP
  - hapmap_to_vcf.py（将hmp.txt转换为vcf）
  - build_plink2_cmd.py（生成plink2的使用指令）
  - origin_maize.geno.selected.hmp.txt（示例基因组文件）
  - origin_maize.geno.selected.vcf（示例转换后的vcf文件）
  - plink2output（用于存放plink2的输出，其内有eigenvec文件）
  
- ### pheno
  - csv_to_tsv.py（将csv转换为tsv）
  - maize.hybrid.train_phe.csv（示例表型组文件）
  - maize.hybrid.train_phe.tsv（示例转换后的tsv文件）

- ### Alignment
  - align_tsv_eigenvec.py（将tsv和eigenvec对齐，即保留相同的个体id并对齐编号）
  - tsv_to_pkl.py（将tsv、csv、eigenvec转化为pkl）
  - maize.hybrid.train_phe.aligned.tsv（对齐后的表型组文件）
  - pca250.aligned.eigenvec（对齐后的基因组文件）
  - pca250.pkl（转化后的基因组文件）

- #### test_output（内为调用DNNGP模型后的输出）

- #### dnngp使用流程图.drawio.png（配合DNNGP官方使用说明更好）
生成模型后，也可使用本项目内的程序来处理基因型文件，使得其可以继续使用DNNGP项目的预测功能





