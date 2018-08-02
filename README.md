##  THIS PAGE UNDER CONSTRUCTION 
Welcome to a tutorial on ClaaTU. This page is under construction to provide a tutorial guide for ClaaTU. The goals of this tutorial is for you to be able to (quickly) understand the theory behind ClaaTU, as well as be able to apply ClaaTU as a tool to your own data. 

### What is CLAATU?
Claatu is a new bioinformatic tool created by Sharpton lab (link here to Chris' website). You can check out our new paper which implements this tool here (add link here). The ClaaTU Algorithm allows us to discover *Cladal Taxonomic Units*, or CTUs. If you've worked doing any microbiome analysis before, you are familiar with an Operational Taxonomic Unit (OTU) table that are output by many softwares such as QIIME, Mothur, or an Amplicon Sequence Variant (ASV) table output by Dada2. ClaaTU can work with any of these common software outputs to discover CTUs in your data. 

### What is a CTU?
A CTU is simply a monophylitic clade of organisms within a phylogentic tree. A CTU can occur at any point along the phylogenetic tree: very close to the tree tips, or very deep within the tree. For the purposes of the tutorial, think of the tree tips as the OTU names. Think of nodes within the tree as CTU names.

[figure here]

### ClaaTU Algorithm
The idea behind ClaaTU is simple. ClaaTU uses two input files:
* An OTU or ASV table 
* A phylogenetic tree representing how OTUs (or ASVs) are related to one another.
* (Optional) a taxonomy file which labels each OTU with taxonomy

The output is a Cladal Taxonomic Unit Matrix m with the following properties:
* the row names correspond to the OTU sample names
* the column names correspond to each node of the phylogenetic tree.
* a m[i,j] corresponds to the sum of counts for each member of clade j found in sample i. 

[figure here]

### Install
First things first, lets install the software. 
```markdown
git git@github.com:arnold3/ClaaTU.git
```

## Dependencies
If you need to, you will have to install the following dependencies
1. [python (2.7.10)](https://www.python.org/downloads/)
2. [dendropy (4.0.2)](https://www.dendropy.org)
3. [scipy (0.13.0b1)](https://www.scipy.org/install.html)
4. [numpy (1.8.0rc1)](https://docs.scipy.org/doc/numpy/user/install.html)


### Files
Lets take a look at what the files are in the ClaaTU folder. Navigate to the location you downloaded Claatu to, then type: 
```markdown
cd Claatu-master
ls
```
You should see the following
```markdown
arc/  bin/  LICENSE.txt  test/  tutorialData/

```
Held within ```bin/``` are all the scripts which will carry out the Claatu algorithm. Held within `tutorialData/` are all of our tutorial files. This is where we are going to be working. Navigate to the `tutorialData` folder 
```markdown 
  cd tutorialData/
  ls
```
Here, we can see there are three files. These are the tree files that you will need to run Claatu. Let's take a look at each one of them in a little more detail.
#### OTU Table

The first file that you will need is the OTU table that is output from QIIME. File format requirements:
- All entries must be tab separated.
- In standard QIIME output, there is text before the OTU table itself. To mark the beginning of the OTU table, QIIME adds a `#OTU ID` before the table. Yes, that is a `#OTU<space>ID<tab>`. This is important to note if you are not getting your OTU table from QIIME, but from another source.
- Column names are *sample IDs*
- Row names are *otu names*

So, now, lets take a look at our particular example:
```markdown
more otu.txt
```
<p align="center">
<kbd>
<img src="https://github.com/arnold3/ClaaTU/blob/master/Claatu-master/tutorialImages/otuTableFormatImage.png" width="50%" height="50%" align="center" style="border:5px solid black"/>
</kbd>
</p>
In our tutorial dataset, we have OTUs sampled from various animals (the column names) and the columns are OTUs. We can see that otu's are named `otu1, otu2, ... otuN`. 

Test image 1
<img src="https://github.com/arnold3/ClaaTU/blob/master/Claatu-master/tutorialImages/otuTableFormatImage.png" width="50%" height="50%" align="center" style="border:5px solid black"/>

#### Bacterial tree
The second required file is a bacterial phylogenetic tree. File format requirements:
- newick format
- The tip labels must match the OTU labels in our otutable (`otu1, otu2, ... otuN`) exactly. Likewise, the names in the `.tre` file must exactly match those in the otu table.

```markdown
more bacterial.tre
```

<p align="center">
<kbd>
<img src="https://github.com/arnold3/ClaaTU/blob/master/Claatu-master/tutorialImages/bacterialTreeFormat.png" width="50%" height="50%" align="center" style="border:5px solid black"/>
</kbd>
</p>

In our tutorial dataset, we see that we have a newick formated file, and that we have tree tips named similarly to the otu table. 

#### OTU Taxonomy table


### Step 1: Prep phylogenetic tree (prep_tree.py)
```python
python ../bin/prep_tree.py bacteria.tre 
```
### Step 2: Get the CTU Matrix (count_tree.py)
### Step 3: Get CTU Stats (clade_stat.py)
### Step 4: Get CTU taxonomy (tax_parser.py)
### Step 5: Node_info.py
### Step 6A: Significance (ptest_tree.py)
### Step 6B: Significance by Group (ptest_tree.py)


You can use the [editor on GitHub](https://github.com/arnold3/CTU_Tutorial/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/arnold3/CTU_Tutorial/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
