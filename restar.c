#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 256

struct Sumatoria
{
	int bloque[MAX];
	int largo;
};
/********************************************************************************************************
	FUNCIONES SECUNDARIAS
********************************************************************************************************/
//Convierte un entero de GMP a una expresion en forma de sumatoria
void convertir(int w, mpz_t numero, struct Sumatoria *S)
{
	mpz_t aux, resto;
	mpz_inits(aux, resto, NULL);
	int i, b, n, largo;	//i = indice; b = 2^W; n = cantidad bloques; largo = cantidad de bits

	mpz_set(aux, numero);				//aux = numero
	largo = mpz_sizeinbase(numero, 2);
	b = (int) pow(2.0, 1.0 * w);
	S->largo = (int) ceil(1.0 * largo / w);	//División techo

	for(i = 0; i < S->largo; i++)
	{
		mpz_fdiv_qr_ui(aux, resto, aux, b);	//aux = aux/b (resto = bloque)
		S->bloque[i] = mpz_get_ui(resto);	//Cada bloque se guarda en una posicion de arreglo
	}
	mpz_clears(aux, resto, NULL);
}

//Inicializa la sumatoria con todos sus valores en cero
void iniciar(struct Sumatoria *A)
{
	for(int i = 0 ; i < MAX ; i++)
		A->bloque[i] = 0;
}

//Impresion en forma de sumatoria
void imprimir(struct Sumatoria A, int b)
{
	int i;
	for(i = 0; i < A.largo ; i++)
	{
		printf("(%d * %d^%d)", A.bloque[i], b, i);
		if(i < A.largo-1) printf(" + ");
	}
	printf("\n");
}

//Busca reducir el largo de una sumatoria para que no existan terminos de la forma (0 * b^i)
void reducir_largo(struct Sumatoria *A)
{
	int i;
	for(i = A->largo - 1 ; i > 0 ; i--)
	{
		if(A->bloque[i] == 0)	A->largo--;
	}
}

/********************************************************************************************************
	OPERACIONES ARITMETICAS
********************************************************************************************************/

//SUMA
void sumar(struct Sumatoria *R, struct Sumatoria s_1, struct Sumatoria s_2, int b)
{
	int i, c, n;
	iniciar(R);	//Inicia el resultado en Cero

	if(s_1.largo > s_2.largo)	n = s_1.largo;	//Determinar largo del mayor sumando para el resultado
	else						n = s_2.largo;

	c = 0;

	for(i = 0; i < n ; i++)	//De 0 a n-1
	{
		R->bloque[i] = (s_1.bloque[i] + s_2.bloque[i] + c) % b;//Z1 <- Xi + Yi + c mod b
		if(s_1.bloque[i] + s_2.bloque[i] + c >= b)	//Xi + Yi + c >= b, entonces
			c = 1;
		else
			c = 0;
		R->bloque[n] =  c;
	}
	if(R->bloque[n] > 0)	R->largo = n + 1;
	else					R->largo = n;
}

//RESTA
void restar(struct Sumatoria *R, struct Sumatoria mi, struct Sumatoria su, int b)
{
	int i, c, n;
	iniciar(R);		//Inicia el resultado en Cero
	n = mi.largo;

	c = 0;

	for(i = 0; i < n ; i++)	//De 0 a n-1
	{
		R->bloque[i] = (mi.bloque[i] - su.bloque[i] + c + 256) % b;	//Z1 <- Xi - Yi + c mod b (+256 = mod para negativos)

		if(mi.bloque[i] - su.bloque[i] + c < 0)	//Xi + Yi + c >= b, entonces
			c = -1;
		else
			c = 0;
	}
	R->largo = n;
	reducir_largo(R);
}

//MULTIPLICACION
void multiplicar(struct Sumatoria *R, struct Sumatoria f_1, struct Sumatoria f_2, int b)
{
	int i, j = 0, c, u, v, uv, m = f_1.largo, n = f_2.largo;
	iniciar(R);	//Inicia el resultado en Cero
	R->largo = n + m;

	for(i = 0 ; i < n ; i++)
	{
		c = 0;
		for(j = 0; j < m ; j++)
		{
			uv = R->bloque[i+j] + f_1.bloque[j] * f_2.bloque[i] + c;	//(uv) =  R{i+j} + f1 * f2 + c
			u = uv / b;
			v = uv % b;
			R->bloque[i+j] = v;
			c = u;
		}
		R->bloque[i+j] = c;				//cambio, mantiene a u aun cuando se cambie de indice
	}
	//R->bloque[n+m-1] = c;	(Segun algoritmo, MALO)
	reducir_largo(R);
}

//CUADRADO
void cuadrado(struct Sumatoria *R, struct Sumatoria S, int b)
{
	int i, j, uv, u, v, c;
	iniciar(R);	//Inicia el resultado en Cero

	for(i = 0 ; i < S.largo ; i++)
	{
		uv = R->bloque[2*i] + S.bloque[i] * S.bloque[i];	// (uv) = R{2i} + Si^2
		u = uv / b;
		v = uv % b;
		R->bloque[2*i] = v;	//R{2i} =v
		c = u;

		for(j = i + 1 ; j < S.largo ; j++)
		{
			uv = R->bloque[i+j] + 2 * S.bloque[i] * S.bloque[j] + c;
			u = uv / b;
			v = uv % b;
			R->bloque[i+j] = v;
			c = u;
		}
		R->bloque[i + S.largo] = u;	//R{2i} =v
	}
	reducir_largo(R);
}

/********************************************************************************************************
	FUNCION PRINCIPAL
********************************************************************************************************/
int main(int argc, char *argv[])
{
	struct Sumatoria A, B, C;
	mpz_t X, Y, R;
	int i, b, n, w;

	mpz_inits(X, Y, R, NULL);
	iniciar(&A);
	iniciar(&B);
	
	mpz_set_str(X, argv[2], 10);
	mpz_set_str(Y, argv[3], 10);
	w = atoi(argv[1]);
	b = (int) pow(2.0, 1.0 * w);

	convertir(w, X, &A);
	convertir(w, Y, &B);

	printf("\n");

	multiplicar(&C, A, B, b);
	imprimir(C, b);

	mpz_pow_ui(R, X, 2);
	gmp_printf("%Zd\n",R);

	mpz_clears(X, Y, R, NULL);
}