## Cose utili
***

Per copiare qualcosa da tolab alla dir corrente:  
  `scp -r <nome_cartella> giovanni.bollotta@tolab.fisica.unimi.it:/path del posto in cui è/` .  
  Attenzione: potrebbe fare casini mettere un argomento di directory prima dei path. Di sicuro questo funziona:   
  `scp -r giovanni.bollotta@tolab.fisica.unimi.it:/path completo che arriva fino alla dire da copiare/ .`

***

Con `quota` visualizzo la quota sul disco
Queste cartelle si possono rimuovere per recuperare spazio su tolab:
`rm -fr .cash .vscode vscodeserver`

***

### Per consegnare:
***
Prima di consegnare la directory dell'esame:   
  clang-format -i *.cpp *.hpp

***
Rimuovere tutto quello che non serve; rinominarla in modo univoco e poi fare  
  `scp -r esame_tnds giovanni.bollotta@tolab.fisica.unimi.it:/home/studenti/giovanni.bollotta/Desktop`
e poi copiarla nel posto giusto da tolab (piger). 
***
### Metodi utili: 
***
Un SelectionSort per `double*`:   

  `void scambiaByReference(double &a, double &b) {
    double c=a;
    a=b;
    b=c;
}

void SelectionSort(double * v, int n){
    for(int k = 0 ; k < n-1 ; k++ ){
        double min = v[k];
        int imin = k;
        for(int j = k+1 ; j < n ; j++ ){
            if(v[j] < min){
                imin = j;
                min = v[j];
            }
        }
        scambiaByReference ( v[k] , v[imin]);
    }
}`  
  
  
Un MergeSort per `double*`:  
  
  `void Merge(double* & v, int left, int right) {
    int mid = left + (right - left)/2;
    int n1 = mid - left + 1;
    int n2 = right - mid;

    double* leftArray = new double[n1];
    double* rightArray = new double[n2];

    for (int i = 0; i < n1; i++)
        leftArray[i] = v[left + i];
    for (int j = 0; j < n2; j++)
        rightArray[j] = v[mid + 1 + j];

    int i = 0, j = 0, k = left;

    // Fusione ordinata
    while (i < n1 && j < n2) {
        if (leftArray[i] <= rightArray[j]) {
            v[k] = leftArray[i];
            i++;
        } else {
            v[k] = rightArray[j];
            j++;
        }
        k++;
    }

    // Copia degli elementi rimanenti
    while (i < n1) {
        v[k] = leftArray[i];
        i++;
        k++;
    }

    while (j < n2) {
        v[k] = rightArray[j];
        j++;
        k++;
    }

    delete[] leftArray;
    delete[] rightArray;
}`
`void MergeSort(double* & v, int left, int right) {
    if (left < right) {
        int mid = left + (right - left)/2;

        MergeSort(v, left, mid);
        MergeSort(v, mid + 1, right);

        Merge(v, left, right);
    }
}`

Si possono modificare abbastanza agevolmente per vector della stl, ovviamente. 
Si possono usare `v.begin()` e `v.end()`, eventualmente modificati con * all'interno della funzione. 





scp -r giovanni.bollotta@tolab.fisica.unimi.it:/home/comune/labTNDS_feb08_compito1 ./Users/Antico/Desktop/Giovanni/TNDS/Temi_esame

