#include <stdio.h>
#include <stdlib.h>

int size_linha=3;
int size_coluna=5;
int minimo;
double *A, *identidade, *Acopia, *matrizFinal;

void estica_encolhe(int indice,double fator){
    for(int i=0; i<minimo; i++)
      identidade[i*minimo + indice]*=(1.0/fator);
       
    for(int i=0; i<size_coluna; i++)
        A[indice*size_coluna + i] *= fator;  
}

void soma_subtracao(int i_static,int j_static, double fator){
    if(j_static<minimo || i_static<minimo){
    for(int i=0; i<minimo; i++)
        identidade[i*minimo + j_static]+=(fator)*identidade[i*minimo + i_static];
    }
    for(int i=0; i<size_coluna; i++)
        A[i_static*size_coluna + i] -= fator* A[j_static*size_coluna + i];  

}

void troca(int i_static, int j_static){
    double aux;
    for(int i=0; i<size_coluna; i++){
        aux=A[i_static*size_coluna + i];
        A[i_static*size_coluna + i]=A[j_static*size_coluna + i];
        A[j_static*size_coluna + i]=aux;
    }

    for(int i=0; i<minimo; i++){
        aux=identidade[i_static + i*minimo];
        identidade[i_static + i*minimo]=identidade[j_static + i*minimo];
        identidade[j_static + i*minimo]=aux;
    }
}

int encontra_pivo(int i,int j){
    int i_troca=-1;
    for(int k = i; k<size_linha; k++){
        if(A[k*size_coluna + j]!=0.0){
            i_troca=k;
            break;
        }
    }
    return i_troca;  
}

void zerar_coluna(int i,int j){  
    for(int k = 0; k<size_linha; k++){
        if (k!=i)
            soma_subtracao(k, i, A[k*size_coluna +j]);
     }
 }

void preencheIdentidade(){
    for(int i=0; i<minimo; i++){
        for(int j=0; j<minimo; j++){
            if(i==j)
                identidade[i*minimo +j]=1.0;
            else
                identidade[i*minimo +j]=0.0;    
            }    
    }

}

void preencheMatrizA(){
    for(int i=0; i<size_linha; i++){
        for(int j=0; j<size_coluna; j++){
            A[i*size_coluna +j]=(j+1*7)/(i+1);//aplicar rand
            Acopia[i*size_coluna +j]=(j+1*7)/(i+1);//aplicar rand                
        }    
    }

}

void  Gauss_Jordan(){
    int i_atual=0;
    int j_atual=0;

    while(i_atual<size_linha && j_atual<size_coluna){
       
        int i_troca = encontra_pivo(i_atual, j_atual);
            
        if(i_troca == -1){
            if(i_atual==size_linha){
                j_atual=size_coluna;
                continue;
            }
            j_atual+=1;
            continue;
        }
      
        if(i_atual!=i_troca)
            troca(i_atual, i_troca);
       
        estica_encolhe(i_atual, 1.0/A[i_atual*size_coluna + j_atual]);
        
        zerar_coluna(i_atual, j_atual);
        
        i_atual+=1;
        j_atual+=1;  
       
    }
   
/* */

 


}

int min(){
    if(size_coluna<size_linha)
        return size_coluna;
    else
         return size_linha;
}

void MultiSequencial(int dim){
 
	float aux=0;
    if(size_coluna<size_linha){
		for(int i=0; i<size_linha; i++){
	      	for(int j=0; j<dim; j++){ 
				matrizFinal[i*dim +j]=0;
                aux=0;
				for(int k=0; k<dim; k++){
					aux+=A[i*dim +k]*identidade[k*dim+j];				
				}
			    matrizFinal[i*dim+j]=aux;
			}
		}
    }
    else{

        for(int i=0; i<minimo; i++){
	      	for(int j=0; j<size_coluna; j++){ 
                matrizFinal[i*size_coluna +j]=0;
                aux=0;
				for(int k=0; k<dim; k++){
					aux+=identidade[i*dim +k]*A[k*size_coluna+j];				
				}
			    matrizFinal[i*size_coluna+j]=aux;
			}
		}

    }
}

int verifica(){
    int aux = 0;
        for(int i=0; i<size_coluna; i++){
            for(int j=0; j<size_linha; j++){
                if(matrizFinal[i*size_coluna+j]!=Acopia[i*size_coluna+j]){
                    aux=1;
                    break;
                }
            }     
            if(aux==1)
            break;            
        }
    return aux;
 }



int main(){
  
    A = (double *) malloc(sizeof(double) * size_coluna*size_linha);
    Acopia  = (double *) malloc(sizeof(double) * size_coluna*size_linha);
    matrizFinal = (double *) malloc(sizeof(double) * size_coluna * size_linha);
    minimo = min();
    identidade= (double *) malloc(sizeof(double) * minimo * minimo);
    preencheIdentidade();
    preencheMatrizA();
    Gauss_Jordan();
    MultiSequencial(min());
    if(verifica(min()))
        printf("\na multiplicação das matrizes obtidas pela fatoração gausseana, é igual a matriz original\\\n");
    free(A);
    free(Acopia);
    free(identidade);
    free(matrizFinal);
    return 0;
}