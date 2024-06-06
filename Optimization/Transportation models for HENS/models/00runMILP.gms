*-------------------------------------------------------------------------------
*--------------   GAMS FILE FOR TRANSPORTATION PROBLEM

*--------------   DEVELOPED BY DR. JUAN M. ZAMORA M.
*--------------   LUIS E. PEDROZA ROBLES A.

*--------------   UNIVERSIDAD AUTÓNOMA METROPOLITANA IZTAPALAPA, MEXICO.
Scalar flag, aux;


$inlinecom /* */
$include options.c

Option MIP = CPLEX;

$include 01input11H-2C.dat


Parameter HRAT, DigHRAT, NO_ITERA;
option seed = 911;
HRAT =       10.000;
DigHRAT =  3;
NO_ITERA =   10;

*----------------------------------------------------------------------------------


*Probabilidad de que un apareamiento sea clonado sobre una corriente madre (0 indica que todos seran clonado, 1 indica que ninguno puede ser clonado)
Scalar ProBC /1/;


$include 02GenIntervalos.dat

Parameter elim;

*-------------------------------------------------------------------------------



Parameter DTML(X,X1), DTMLH(X1), DTMLC(X);


Loop[(X,X1)$(Ord(X) LE (NOK) and Ord(X1) LE (NOK) and Ord(X1) GE Ord(X)),
  If((TH(X)-TC(X1)) EQ (TH(X+1)-TC(X1+1)),
    DTML(X,X1) = ((TH(X)-TC(X1))+(TH(X+1)-TC(X1+1)))/2;
  else 
    DTML(X,X1) = ((TH(X)-TC(X1))-(TH(X+1)-TC(X1+1)))/LOG((TH(X)-TC(X1))/(TH(X+1)-TC(X1+1)));
  );
];
DTML(K,K1)$(OrdK(K) GT OrdK1(K1)) = 1e-15;


Loop[X$(Ord(X) LE (NOK) and TH(X) GT TTCU and TH(X+1) GT TSCU),
  If((TH(X)-TTCU) EQ (TH(X+1)-TSCU),
    DTMLC(X) = ((TH(X)-TTCU)+(TH(X+1)-TSCU))/2;
  else 
    DTMLC(X) = ((TH(X)-TTCU)-(TH(X+1)-TSCU))/LOG((TH(X)-TTCU)/(TH(X+1)-TSCU));
  );
];
DTMLC(K)$(DTMLC(K) EQ 0) = 1e-15;


Loop[X1$(Ord(X1) LE (NOK) and TSHU GT TC(X1) and TTHU GT TC(X1+1)),
  If(TSHU-(TC(X1)) EQ (TTHU-TC(X1+1)),
    DTMLH(X1) = (TSHU-(TC(X1))+(TTHU-TC(X1+1)))/2;
  else 
    DTMLH(X1) = ((TSHU-TC(X1))-(TTHU-TC(X1+1)))/LOG((TSHU-(TC(X1)))/(TTHU-TC(X1+1)));
  );
];
DTMLH(K)$(DTMLH(K) EQ 0) = 1e-15;


Parameter U(I,J), Uh(J), Uc(I);
U(I,J) = 1/(1/FHTCI(I)+1/FHTCJ(J));
Uh(J) = 1/(1/FHTCHU+1/FHTCJ(J));
Uc(I) = 1/(1/FHTCI(I)+1/FHTCCU);

Display DTML, DTMLC, DTMLH, U, Uh, Uc;

Parameter M(I,J), MHU(J), MCU(I);
M(I,J) = min(sum(K,HS(I,K)),sum(K1,HD(J,K1)));
MHU(J) = sum(K1,HD(J,K1));
MCU(I) = sum(K,HS(I,K));


*-------------------------------------------------------------------------------
*--------------  MODELO DE TRANSPORTE
Parameter qhut Servicios totales de calentamiento
          qcut Servicios totales de enfriamiento
          qrec Calor total recuperado ;

Variable funobj;
Positive variables
           q(I,X,J,X1) Calor recuperado
           qhu(J,X1)   Calor proveniente de servicios de calentamiento
           qcu(I,X)    Calor proveniente de servicios de enfriamiento;

Binary variables Y(I,J), YHU(J), YCU(I);           

Equation obj función a minimizar;
obj..

funobj =E= Sum((I,J),CCFXCGL*Y(I,J))+Sum((I,K,J,K1),CCPRECL*(q(I,K,J,K1)/(U(I,J)*DTML(K,K1))))+
           Sum(J,CCFXCGHUL*YHU(J))+Sum((J,K1),CCPRECHUL*(qhu(J,K1)/(Uh(J)*DTMLH(K1))))+
           Sum(I,CCFXCGCUL*YCU(I))+Sum((I,K),CCPRECCUL*(qcu(I,K)/(Uc(I)*DTMLC(K))))+
           Sum((J,K1),HUCOST*qhu(J,K1))+Sum((I,K),CUCOST*qcu(I,K));
 

Equation SUPPLY(I,X) Restricciones de suministro;
SUPPLY(I,K).. Sum((J,K1)$(OrdK1(K1) GE OrdK(K)),q(I,K,J,K1))+qcu(I,K) =E= HS(I,K);

Equation DEMAND(J,X1) Restricciones de demanda;
DEMAND(J,K1).. Sum((I,K)$(OrdK(K) LE OrdK1(K1)),q(I,K,J,K1))+qhu(J,K1) =E= HD(J,K1);

Equation CotasQ(I,J) Cotas para la transferencia de calor entre corrientes de proceso;
CotasQ(I,J).. Sum((K,K1)$(OrdK1(K1) GE OrdK(K)) , q(I,K,J,K1)) =L= M(I,J)*Y(I,J);

Equation CotasQHU(J) Cotas para la transferencia de calor de corrientes frias con servicios de calentamiento;
CotasQHU(J).. Sum(K1, qhu(J,K1)) =L= MHU(J)*YHU(J);

Equation CotasQCU(I) Cotas para la transferencia de calor de corrientes calientes con servicios de enfriamiento;
CotasQCU(I).. Sum(K, qcu(I,K)) =L= MCU(I)*YCU(I);






Model transport /all/;

transport.prioropt = 0;   /*Use priorities. Priorities should be assigned based on your knowledge of the problem*/

transport.iterlim = 1e9;    /*Sets the simplex iteration limit. Simplex algorithms will terminate and pass on the current solution to GAMS.*/
transport.reslim = 10000;      /*Sets the time limit in seconds. The algorithm will terminate and pass on the current solution to GAMS.*/
transport.workspace = 10000; /*Request of memory*/

transport.optca = 1e-7;   /*Absolute optimality criterion for a MIP problem. The OptCA option asks Cplex to stop when ABS(ObjFun-BestPossibleIntegerSol) < OptCa*/
transport.optcr = 1e-9;   /*Relative optimality criterion for a MIP problem. The OptCR option asks Cplex to stop when |BP-BF|/(qe-10+|BF|)<OptCR*/


Solve transport using MIP minimizing funobj;

