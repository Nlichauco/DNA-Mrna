
#Summary: This Program is meant to take the same gene sequence and put it through
#three reading frames in order to figure out which is best for reading. 

#The input is a gene, the gene then goes through the reading frames and is transcribed
#and then translated into protein, all within the program.

#The output: In each reading frame the gene is changed and broken up into codons,
#each codon is matched to the according protein and that is printed out.
#if there are any starts or stops inbetween it will be pointed out in the specific reading frame.

#The intended purpose is to give all options for reading the gene and let the 
#user decide which is most effective for their purpose. 

import BioDNA

def main():
    
    print ("\n       ---+++ Central Dogma +++---\n")
    
    AAtable = BioDNA.makeAminoAcidTable()
    
    for (tripleRNA, AAcode) in AAtable.items():
        print (tripleRNA, AAcode)

        
    DNA = "ATGTGTACGCAAATGATATCGTATTAG";

    print("-"*40)
    print (" DNA: 5'", DNA, "3'")
    
    # --------- TRANSCRIBE --------------------------
    
    mRNA=BioDNA.transcribe(DNA)
    print("        ", "|"*len(mRNA))
    print("mRNA: 5'", mRNA, "3'")
  

    print("\n==============================================================")
    print("Reading Frame #1")
    # --------- TRANSLATE Reading Frame #1 --------------------------
    p=BioDNA.translate(mRNA, AAtable)
    print("RNA length:",len(mRNA))
    print("RNA:    ",mRNA,)
    print("        "," | "*(len(mRNA)//3))
    print("Protein:", p,"\n")
    #most blocks of the if statements are the same throughout the readings
    #the only thing that changes is the gene.
    stop=p.find("***")
    if stop == len(p)-3:
        print("There is a stop codon at the end")
    
    else:
        print("There is no stop codon")
        
    #this searches for the index of the start
    find=mRNA.find("AUG")
    if find == 0:
        print("There is a met at", find)
    
    else:
        print("There is no Met at the start")
    #below I look through the middle 
    
    midend=len(p)-3
    mid=p[3:midend]
    mid_start=mid.find("Met")
    if mid_start >=0:
        print("there is a Met inbetween at position,", mid_start+3)
    else:
        print("there is no Met inbetween")
        
    mid_end=mid.find("***")
    if mid_end >=0:
        print("there is a *** inbetween at position,", mid_end+3)
    else:
        print("there is no *** inbetween")
    
    
    print("-"*40)   
    
    print("==============================================================\n")
    
    
    print("\n==============================================================")
    print("Reading Frame #2")
    # --------- TRANSLATE Reading Frame #2 --------------------------
    #Essentially the same codes as the last one but specific to this gene
    mRNA2=mRNA[1:]
    e=BioDNA.translate(mRNA2, AAtable)
    print("RNA length: ", len(mRNA2),"bp")
    print("RNA:  5'", mRNA2, "3'")
    print("        "," | "*(len(mRNA2)//3))
    print("Protein:",e,"\n")
    find2=e.find("Met")
    end=e.find("***")
    if end == len(p)-3:
        print("There is a stop at the end")
    else:
        print("There is no stop at the end") 
        
    
    if find2 == 0:
        print("There is a Met at", find2)
    
    else:
        print("There is no Met at the start")
    #ALl the text below till the next reading frame is specifically for the middle
    #of the code.
    midend1=len(e)-3
    mid1=e[3:midend1]
    mid_start1=mid1.find("Met")
    #Checking the middle below
    
    if mid_start1 >=0:
        print("there is a start inbetween at position,", mid_start1+3)
    else:
        print("there is no Met inbetween")
    mid_end1=mid1.find("***")
    if mid_end1 >=0:
        print("there is a *** inbetween at position,", mid_end1+3)
    else:
        print("there is no *** inbetween")
    print("==============================================================\n")
    
    
    print ("\n==============================================================")
    print ("Reading Frame #3")
    # --------- TRANSLATE Reading Frame #3 --------------------------
    mRNA3=mRNA2[1:]
    len3=len(mRNA3)
    print("RNA length: ",len3,"bp")
    
    t=BioDNA.translate(mRNA3, AAtable)
    mRNA3=BioDNA.transcribe(mRNA3)
    print("RNA:    ",mRNA3)
    print("        ", " | "*(len(mRNA3)//3))
    print("Protein:",t,"\n")
    find3=t.find("Met")
    finish=t.find("***")
    if finish == len(t)-3:
        print("There is a stop at the end")
    else:
        print("There is no stop codon")
        
    
    if find3 == 0:
        print("There is a met at", find3)
    else:
        print("There is no Met at the start")
        
    midend2=len(t)-3
    mid2=t[3:midend2]
    mid_start2=mid2.find("Met")
    if mid_start2 >=0:
        print("there is a Met inbetween at position,", mid_start2+3)
    else:
        print("There is no Met inbetween")
    mid_end2=mid2.find("***")
    if mid_end2 >=0:
        print("there is a *** inbetween at position,", mid_end2+3)
    else:
        print("there is no *** in the middle")
    
                 
    
    
    
    
    print ("==============================================================\n")
    
    

        
# --- end main() ---------


#---------------------------------------------------------
# Python starts here ("call" the main() function at start
if __name__ == '__main__':
    main()
#---------------------------------------------------------  
    
