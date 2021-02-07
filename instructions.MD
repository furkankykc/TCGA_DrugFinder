2) You will retrieve top 10 corresponding PubChem CIDs(Compound IDs) for each drug name using 
   * eutils: https://www.ncbi.nlm.nih.gov/pccompound/?term=%22aspirin%22%5BSynonym%5D
* PubChem Compounds might not have much readable description;
   * therefore you will efetch the MeSH term of each compound CID. 
   * Combine all info into a 6 dimensional hash of (Patient barcode, CID, immunotherapy/chemotherapy,smoking/nonsmoking,Mesh Term for the CID, survival days). 
   * For each drug Mesh Term, create a summary output of how many patients survived, how many of them were smoking, immunotherapy/chemotherapy/targeted molecular therapy in decreasing order of mean survival  
   * (you may use multi-dimensional hash and sort by value of the hash(survival). 
   * Your program should then ask for user input of two compound names and print their statistics side by side* Write a program to run this task both from commandline and through a user-interface.  