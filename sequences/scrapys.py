import json
from Bio import Entrez, SeqIO
from subprocess import run, CalledProcessError


Entrez.email = "joserodolfo.silva@ufpe.br"
# Comando para baixar sequências de nucleotídeos
def getsequences(db: str, id: str, rettype: str, retmode:str ):
    handle = Entrez.efetch(db=db,id=id,rettype=rettype,retmode=retmode)
    with open("NZ_CP041838.1.fasta","w") as arquivo:
        arquivo.write(handle.read())
    handle.close()


def get_read_id(term:str, #Especificar a pesquisa que está sendo feita
                retmax:int,
                db='sra'):
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

    
# Comanando para baixar informações do arquivo SRA  
def get_reads(ids:str, db='sra'):
    handle = Entrez.efetch(db=db, id=ids, rettype='fasta', retmode='text')
    sequences = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return sequences

# Comando usado para baixar informações de variantes genéticas
def get_genvariants_info(db='snp',term=str, retmax=int):
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    sequencias = Entrez.read(handle)
    handle.close()
    return sequencias


# Baixando variantes genéticas
def get_genvariants(db='snp', id=int, rettype="gb", retmode="text"):
    handle = Entrez.efetch(db=db, id=id, rettype="gb", retmode=retmode)

# Comando para baixar as reads
def get_reads_cli(serial_number):
    try:
        command = ['prefetch', serial_number]
        
        result = run(command, capture_output=True, text=True, check=True)
        
        print("Prefetch executado com sucesso!")
        print("Saída:", result.stdout)
    
    except CalledProcessError as e:
        print("Erro ao executar prefetch:")
        print(e.stderr)


