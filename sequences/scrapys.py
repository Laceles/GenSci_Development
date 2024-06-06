import json
from Bio import Entrez

def getseq(db: str, id: str, rettype: str, retmode:str ):
    Entrez.email = "joserodolfo.silva@ufpe.br"
    handle = Entrez.efetch(db=db,id=id,rettype=rettype,retmode=retmode)
    with open("NZ_CP041838.1.fasta","w") as arquivo:
        arquivo.write(handle.read())
    handle.close()