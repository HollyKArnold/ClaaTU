##  THIS PAGE UNDER CONSTRUCTION 
Welcome to a tutorial on ClaaTU. This page is under construction to provide a tutorial guide for ClaaTU. The goal of this tutorial is for you to be able to understand the theory behind ClaaTU, as well as be able to apply ClaaTU as a tool to your own data. 

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
git 
```

### Step 1: Prep phylogenetic tree (prep_tree.py)
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
