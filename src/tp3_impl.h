#ifndef TP3_IMPL_H
#define TP3_IMPL_H

#include "tp3.h"

#include <limits>
#include <algorithm>
#include <iostream>
#include <map>

///////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 1
////
template <typename iterador>
inline void obtener_size_distancia(int& size, int& distancia, iterador input_begin,
                                   iterador input_end){
    int maximo = int(*input_begin);
    int minimo = int(*input_begin);
    for (auto it = input_begin; it != input_end; ++it) {
        if (maximo < int(*it))
            maximo = int(*it);
        else if (int(*it) < minimo)
            minimo = int(*it);
    }

    size = maximo - minimo + 1;
    distancia = 0 - minimo;
}

template <typename iterador, typename bucket>
vector<bucket> generar_buckets(iterador input_begin, iterador input_end) {
    if(input_begin != input_end) {
        int size, distancia;
        //Obtengo el tamanio del vector de buckets
        obtener_size_distancia(size, distancia, input_begin, input_end);
        vector<bucket> buckets;
        buckets.resize(size);

        //Inserto elementos en sus buckets correspondientes
        for (auto it = input_begin; it != input_end; ++it) {
            bucket *b = &buckets[int(*it) + distancia];
            (*b).insert((*b).end(), *it);
        }

        return buckets;

    } else return vector<bucket>();
}

template <typename bucket>
vector<typename bucket::value_type> aplanar_buckets(const std::vector<bucket> & B) {
    vector<typename bucket::value_type> res;
    for(auto b: B)
        for(auto c: b) {
            res.push_back(c);
        }
    return res;
}

///////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 2
////

fajo ordenar_por_probabilidad(const fajo& falsos_conocidos, const fajo & a_ordenar) {
    // conteo de cantidad de falsos por año
    std::map<int, std::set<billete> > falsos_por_anio;
    for (auto & b: falsos_conocidos)    //Itera F veces
        falsos_por_anio[int(b)].insert(b); //O( log(M) )
    //=> O( F*log(M) )

    fajo resultado;
    resultado.reserve(a_ordenar.size());
    // clasifico billetes a ordenar
    for (auto & b : a_ordenar) //Itera N veces
        if (falsos_por_anio[int(b)].count(b)) {     //O( log(#anios) + log(M) ) = O( log(M) ) porque los anios estan acotados por M
            resultado.push_back(b); //O(1)
        } else {
            resultado.push_back(billete(b.numero_de_serie, falsos_por_anio[int(b)].size())); //O(1)
        }
    //=> O( N*log(M) )

    std::sort(resultado.rbegin(), resultado.rend());//O(N log N)

    //=> O( F*log(M) ) + O( N*log(M) ) + O( N*log(N) ) = O( F*log(M) + N*(log(M)+log(N)) )
    // = O( F*log(M) + N*log(M*N) )

    return resultado;
}

///////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 3
////

inline Matriz Bloque11(const Matriz& M){
    unsigned long size = M.size()/2;
    Matriz res(size, vector<double>(size,0.0));
    for(unsigned long i=0; i< size; i++){
        for(unsigned long j=0; j< size; j++){
            res[i][j] = M[i][j];
        }
    }
    return res;
}

inline Matriz Bloque12(const Matriz& M){
    unsigned long size = M.size()/2;
    Matriz res(size, vector<double>(size,0.0));
    for(unsigned long i=0; i< size; i++){
        for(unsigned long j=0; j< size; j++){
            res[i][j] = M[i][j+size];
        }
    }
    return res;
}

inline Matriz Bloque21(const Matriz& M){
    unsigned long size = M.size()/2;
    Matriz res(size, vector<double>(size,0.0));
    for(unsigned long i=0; i< size; i++){
        for(unsigned long j=0; j< size; j++){
            res[i][j] = M[i+size][j];
        }
    }
    return res;
}

inline Matriz Bloque22(const Matriz& M){
    unsigned long size = M.size()/2;
    Matriz res(size, vector<double>(size,0.0));
    for(unsigned long i=0; i< size; i++){
        for(unsigned long j=0; j< size; j++){
            res[i][j] = M[i+size][j+size];
        }
    }
    return res;
}

inline Matriz SumarMatrices(const Matriz& A, const Matriz& B){
    Matriz res(A.size(), vector<double>(A.size(),0.0));
    for(unsigned long i=0; i< A.size(); i++){
        for(unsigned long j=0; j< A.size(); j++){
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

inline Matriz RestarMatrices(const Matriz& A, const Matriz& B){
    Matriz res(A.size(), vector<double>(A.size(),0.0));
    for(unsigned long i=0; i< A.size(); i++){
        for(unsigned long j=0; j< A.size(); j++){
            res[i][j] = A[i][j] - B[i][j];
        }
    }
    return res;
}

inline Matriz ArmarMatrizDesdeBloques(const Matriz& C11, const Matriz& C12, const Matriz& C21,
                                      const Matriz& C22){
    unsigned long size = C11.size()*2;
    Matriz res(size, vector<double>(size,0.0));
    for(unsigned long i=0; i< size/2; i++){
        for(unsigned long j=0; j< size/2; j++){
            res[i][j] = C11[i][j];
            res[i][j+C11.size()] = C12[i][j];
            res[i+C11.size()][j] = C21[i][j];
            res[i+C11.size()][j+C11.size()] = C22[i][j];
        }
    }
    return res;
}

inline Matriz multiplicar_strassen(const Matriz& A, const Matriz& B, int K) {
    if(int(A.size()) <= K) return multiplicar(A,B);

    /** SEPARO AMBAS MATRICES EN BLOQUES */
    Matriz A11 = Bloque11(A);
    Matriz A12 = Bloque12(A);
    Matriz A21 = Bloque21(A);
    Matriz A22 = Bloque22(A);
    Matriz B11 = Bloque11(B);
    Matriz B12 = Bloque12(B);
    Matriz B21 = Bloque21(B);
    Matriz B22 = Bloque22(B);

    /** ARMO LAS 7 MATRICES */
    Matriz M1_A = SumarMatrices(A11, A22);
    Matriz M1_B = SumarMatrices(B11, B22);
    Matriz M2_A = SumarMatrices(A21, A22);
    Matriz M2_B = B11;
    Matriz M3_A = A11;
    Matriz M3_B = RestarMatrices(B12, B22);
    Matriz M4_A = A22;
    Matriz M4_B = RestarMatrices(B21, B11);
    Matriz M5_A = SumarMatrices(A11, A12);
    Matriz M5_B = B22;
    Matriz M6_A = RestarMatrices(A21, A11);
    Matriz M6_B = SumarMatrices(B11, B12);
    Matriz M7_A = RestarMatrices(A12, A22);
    Matriz M7_B = SumarMatrices(B21, B22);

    Matriz M1 = multiplicar_strassen(M1_A, M1_B, K);
    Matriz M2 = multiplicar_strassen(M2_A, M2_B, K);
    Matriz M3 = multiplicar_strassen(M3_A, M3_B, K);
    Matriz M4 = multiplicar_strassen(M4_A, M4_B, K);
    Matriz M5 = multiplicar_strassen(M5_A, M5_B, K);
    Matriz M6 = multiplicar_strassen(M6_A, M6_B, K);
    Matriz M7 = multiplicar_strassen(M7_A, M7_B, K);

    /** ARMO LOS BLOQUES DE LA MATRIZ RESULTADO */
    Matriz C11 = SumarMatrices(M1, M4);
    C11 = RestarMatrices(C11, M5);
    C11 = SumarMatrices(C11, M7);
    Matriz C12 = SumarMatrices(M3, M5);
    Matriz C21 = SumarMatrices(M2, M4);
    Matriz C22 = RestarMatrices(M1, M2);
    C22 = SumarMatrices(C22, M3);
    C22 = SumarMatrices(C22, M6);

    /** JUNTO LOS BLOQUES RESULTANTES EN UNA UNICA MATRIZ */
    return ArmarMatrizDesdeBloques(C11, C12, C21, C22);
}

#endif // TP3_IMPL_H
