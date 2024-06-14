import json
from Bio import Entrez, SeqIO


user = Entrez.email = "joserodolfo.silva@ufpe.br"
# Comando para baixar sequências de nucleotídeos
def getsequences(db: str, id: str, rettype: str, retmode:str ):
    handle = Entrez.efetch(db=db,id=id,rettype=rettype,retmode=retmode)
    with open("NZ_CP041838.1.fasta","w") as arquivo:
        arquivo.write(handle.read())
    handle.close()
    
# Comanando para baixar informações do arquivo SRA  
def get_reads(ids, db='nucleotide'):
    handle = Entrez.efetch(db=db, id=ids, rettyoe='fasta', retmode='text')
    sequences = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return sequences