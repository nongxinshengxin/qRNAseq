# qRNAseq
<a name="kz53z"></a>
## 写在前面
前期我们开发了可以便捷完成RNAseq下游分析的界面化shiny APP RNAdiff，该APP可以对已经完成比对以及表达矩阵计算的数据进行差异表达基因分析、GO/KEGG富集分析、TPM计算、火山图热图气泡图绘制等【】，但很多时候，公司返回的下机数据往往只是经过简单质控的fastq文件，如何从fastq文件得到表达矩阵，来方便我们使用RNAdiff进行进一步分析呢？这里就涉及到了RNAseq上游分析流程，我们在获得下机数据后（一般公司会提供质控后的数据），需要用比对软件将数据比对到参考基因组上，获得bam文件，再进一步使用定量工具软件获得表达量矩阵。**hisat2+featureCount**是常见的上游分析工具，hisat2为比对软件，featureCounts为定量工具，关于这两个软件的使用教程网络上非常全面，但面对大量样本的转录组数据时，重复调用这些软件会非常麻烦且消耗时间。今天我们将为大家分享一个**超便捷进行RNAseq上游分析的pipeline脚本** qRNAseq.py（_quick RNAseq analysis pipeline_）。只需**一行代码**，即可完成从下机数据到表达矩阵的全套RNAseq上游分析。

<a name="YktoE"></a>
## 环境配置
使用该脚本前，需要安装hisat2、featureCounts和samtools。教程很多，这里建议用conda安装。
```shell
#可以先创建一个专门进行RNAseq分析的环境
conda creat -n rnaseq
#激活环境
conda activate rnaseq
#安装hisat2
conda install hisat2
#安装samtools
conda install samtools
#安装featureCounts，featureCounts被包装在subread软件中
conda install subread
```
<a name="e6ZZS"></a>
## 获取源码
脚本已上传GitHub仓库，可以通过以下链接获取

<a name="fTmQP"></a>
## 如何使用
在使用本脚本前，建议将脚本文件添加至环境变量，如未添加至环境变量，使用脚本时需保证脚本在当前目录下或直接调用脚本的绝对路径。下面所有示例代码，脚本均已添加至环境变量。<br />本文测试所用软件版本为hisat2 v2.2.1，samtools v1.10，featureCounts v2.0.1。

<a name="Ll6k3"></a>
### 使用说明
脚本的使用说明如下：
```markdown
usage: qRNAseq.py [-h] [-d DIRECTORY] [-g GENOME] [-a ANNOTATION] [-o OUTPUT] [-r READ] [-s STRAND] [--intro INTRO] [--thread THREAD]

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        directory that containing FASTQ file(s), .fq/.fastq(.gz) fq文件所在文件夹
  -g GENOME, --genome GENOME
                        reference genome 参考基因组路径及参考基因组索引文件的前缀
  -a ANNOTATION, --annotation ANNOTATION
                        annotation file in GTF format 注释GTF路径
  -o OUTPUT, --output OUTPUT
                        output file path 结果输出文件存放路径
  -r READ, --read READ  Single-read or Paired-End 单端测序或双端测序 -r p Paired-End [default] 双端测序(默认) -r s Single-read 单端测序
  -s STRAND, --strand STRAND
                        strand-specific information 链特异性信息 -s 1 stranded [default] 链特异性(默认) -s 0 unstranded 非链特异性
  -i INTRO, --intro INTRO
                        max intron length, default 4000 bp 最大内含子长度，默认4000bp
  -t THREAD, --thread THREAD
                        Running threads, default 1 运行线程，默认1
```
使用如下命令可以查看帮助：
```markdown
qRNAseq.py -h [--help]
```
该脚本有四个必需参数，分别是<br />-d (--directory) 后接所有进行分析的下机数据所在目录的路径。<br />-g (--genome) 后接参考基因组所在路径，或者是根据参考基因组构建的hisat2索引所在路径及其前缀，如果你提供的路径仅仅是参考基因组文件路径，没有索引，脚本会自动构建索引文件，可以正常后续运行。<br />-a (--annotation) 后接基因结构注释的GTF文件所在路径。注意，一定是GTF格式文件，因为featureCounts只接受GTF文件。<br />-o (--output) 后接最终输出结果存放的目录。目录下会产生两个新文件夹，bam文件夹用来存储比对后的bam文件，featureCounts文件夹用来存储表达矩阵。

<a name="zoLBt"></a>
### 命令行界面
基础命令：
```shell
qRNAseq.py -d /your/path/to/fqdirectory -g /your/path/to/genome.fasta -a /your/path/to/reference.GTF -o /your/path/to/outputdirectory
```

-r (--read)参数决定输入的数据是单端测序(Single-read)还是双端测序(Paired-End)，默认是双端测序。如你的数据为单端测序数据，命令行如下：
```shell
qRNAseq.py -d /your/path/to/fqdirectory -g /your/path/to/genome.fasta -a /your/path/to/reference.GTF -o /your/path/to/outputdirectory -r s
```

-s (--strand)参数决定输入的数据是链特异性RNAseq还是非链特异性RNAseq，默认是链特异性RNAseq。如你的数据是非链特异性RNAseq，命令行如下：
```shell
qRNAseq.py -d /your/path/to/fqdirectory -g /your/path/to/genome.fasta -a /your/path/to/reference.GTF -o /your/path/to/outputdirectory -s 0
```

-i (--intro)参数控制hisat2软件中最大内含子长度，默认是4000bp。修改该参数命令行如下：
```shell
qRNAseq.py -d /your/path/to/fqdirectory -g /your/path/to/genome.fasta -a /your/path/to/reference.GTF -o /your/path/to/outputdirectory -i 5000
```

-t (--thread)参数决定软件运行的线程数，默认一个线程。修改该参数命令行如下：
```shell
qRNAseq.py -d /your/path/to/fqdirectory -g /your/path/to/genome.fasta -a /your/path/to/reference.GTF -o /your/path/to/outputdirectory -t 16
```

<a name="QzU1T"></a>
### 输出结果
最终输出的结果为比对后的bam文件和定量后的表达矩阵。
输出的bam文件可以IGV可视化，输出的表达矩阵可以使用我们开发的RNAdiff APP进行后续下游分析。

<a name="TsREY"></a>
## 脚本价值
qRNAseq 可以说是RNAdiff的BEST Match，包括了RNAseq分析的全套流程，且极大降低了学习和使用门槛，如果大家能很好的结合使用qRNAseq 和RNAdiff，将会真正实现RNAdseq分析自由，再也不需要在公司花费以样本数为单位计算的高昂分析费用。**毕竟，不管什么时候，将命运掌握在自己手里，才是最让人安心的。**



