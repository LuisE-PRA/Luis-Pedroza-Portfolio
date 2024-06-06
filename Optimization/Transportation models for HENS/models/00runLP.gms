
*-------------------------------------------------------------------------------
*--------------   GAMS FILE FOR TRANSPORTATION PROBLEM

*--------------   DEVELOPED BY DR. JUAN M. ZAMORA M.
*--------------   LUIS E. PEDROZA ROBLES A.

*--------------   UNIVERSIDAD AUTÓNOMA METROPOLITANA IZTAPALAPA, MEXICO.
Scalar flag, aux;


$inlinecom /* */
$include options.c

Option LP = CPLEX;

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
DTMLH(K1)$(DTMLH(K1) EQ 0) = 1e-15;


Parameter U(I,J), Uh(J), Uc(I);
U(I,J) = 1/(1/FHTCI(I)+1/FHTCJ(J));
Uh(J) = 1/(1/FHTCHU+1/FHTCJ(J));
Uc(I) = 1/(1/FHTCI(I)+1/FHTCCU);


Parameter M(I,J), MHU(J), MCU(I);
M(I,J) = min(sum(K,HS(I,K)),sum(K1,HD(J,K1)));
MHU(J) = sum(K1,HD(J,K1));
MCU(I) = sum(K,HS(I,K));

Scalar Ep /1e-6/;
Ep = Ep+1;

Parameters AL(I,X,J,X1), CL(I,X,J,X1), CL2(I,X,J,X1), mL(I,X,J,X1), bL(I,X,J,X1), ML(I,X,J,X1);
Parameters ALHU(J,X1), CLHU(J,X1), CL2HU(J,X1), mLHU(J,X1), bLHU(J,X1), MLHU(J,X1);
Parameters ALCU(I,X), CLCU(I,X), CL2CU(I,X), mLCU(I,X), bLCU(I,X), MLCU(I,X);

Loop[(i,k,j,k1)$(HS(I,K) GT 0 and HD(J,K1) GT 0 and DTML(K,K1) GT 1e-10),
    AL(I,K,J,K1) = min(HS(I,K),HD(J,K1))/(U(I,J)*DTML(K,K1));
    CL(I,K,J,K1) = FANNUAL*(CCFXCG+CCPREC*AL(I,K,J,K1)**CCAEXP);
    CL2(I,K,J,K1) = FANNUAL*(CCFXCG+CCPREC*(AL(I,K,J,K1)/Ep)**CCAEXP);
    mL(I,K,J,K1) = (CL(I,K,J,K1)-CL2(I,K,J,K1))/(AL(I,K,J,K1)-(AL(I,K,J,K1)/Ep));
    bL(I,k,J,k1) = -mL(I,K,J,K1)*AL(I,K,J,K1)+CL(I,K,J,K1);
];

Loop[(j,k1)$(HD(J,K1) GT 0 and DTMLH(K1) GT 1e-10),
    ALHU(J,K1) = HD(J,K1)/(Uh(J)*DTMLH(K1));
    CLHU(J,K1) = FANNUAL*(CCFXCGHU+CCPRECHU*ALHU(J,K1)**CCAEXPHU);
    CL2HU(J,K1) = FANNUAL*(CCFXCGHU+CCPRECHU*(ALHU(J,K1)/Ep)**CCAEXPHU);
    mLHU(J,K1) = (CLHU(J,K1)-CL2HU(J,K1))/(ALHU(J,K1)-(ALHU(J,K1)/Ep));
    bLHU(J,k1) = -mLHU(J,K1)*ALHU(J,K1)+CLHU(J,K1);
];

Loop[(i,k)$(HS(I,K) GT 0 and DTMLC(K) GT 1e-10),
    ALCU(I,K) = HS(I,K)/(Uc(I)*DTMLC(K));
    CLCU(I,K) = FANNUAL*(CCFXCGCU+CCPRECCU*ALCU(I,K)**CCAEXPCU);
    CL2CU(I,K) = FANNUAL*(CCFXCGCU+CCPRECCU*(ALCU(I,K)/Ep)**CCAEXPCU);
    mLCU(I,K) = (CLCU(I,K)-CL2CU(I,K))/(ALCU(I,K)-(ALCU(I,K)/Ep));
    bLCU(I,k) = -mLCU(I,K)*ALCU(I,K)+CLCU(I,K);
];

*-------------------------------------------------------------------------------
*--------------  MODELO DE TRANSPORTE MILP
Parameter qhut Servicios totales de calentamiento
          qcut Servicios totales de enfriamiento
          qrec Calor total recuperado ;

Variable funobj, funobjL;
Positive variables
           q(I,X,J,X1) Calor recuperado
           qhu(J,X1)   Calor proveniente de servicios de calentamiento
           qcu(I,X)    Calor proveniente de servicios de enfriamiento
           YL(I,X,J,X1), YHUL(J,X1), YCUL(I,X);

Binary variables Y(I,J), YHU(J), YCU(I);           
          

*-------------------------------------------------------------------------------
*--------------  MODELO DE TRANSPORTE LP

Equation objT función a minimizar;
objT..

funobjL =E=  Sum((J,K1),HUCOST*qhu(J,K1))+Sum((I,K),CUCOST*qcu(I,K));

*-------------------------------------------------------------------------------


Equation obj función a minimizar;
obj..
funobj =E= Sum((I,K,J,K1),bL(I,K,J,K1)*YL(I,K,J,K1))+Sum((I,K,J,K1),mL(I,K,J,K1)*(q(I,K,J,K1)/(U(I,J)*DTML(K,K1))))+
           Sum((J,K1),bLHU(J,K1)*YHUL(J,K1))+Sum((J,K1),mLHU(J,K1)*(qhu(J,K1)/(Uh(J)*DTMLH(K1))))+
           Sum((I,K),bLCU(I,K)*YCUL(I,K))+Sum((I,K),mLCU(I,K)*(qcu(I,K)/(Uc(I)*DTMLC(K))))+
           Sum((J,K1),HUCOST*qhu(J,K1))+Sum((I,K),CUCOST*qcu(I,K));

 

Equation SUPPLY(I,X) Restricciones de suministro;
SUPPLY(I,K).. Sum((J,K1)$(OrdK1(K1) GE OrdK(K)),q(I,K,J,K1))+qcu(I,K) =E= HS(I,K);

Equation DEMAND(J,X1) Restricciones de demanda;
DEMAND(J,K1).. Sum((I,K)$(OrdK(K) LE OrdK1(K1)),q(I,K,J,K1))+qhu(J,K1) =E= HD(J,K1);


Equation CotasQ(I,X,J,X1) Cotas para la transferencia de calor entre corrientes de proceso;
CotasQ(I,K,J,K1)..  q(I,K,J,K1) =L= M(I,J)*YL(I,K,J,K1);

Equation CotasQHU(J,X1) Cotas para la transferencia de calor de corrientes frias con servicios de calentamiento;
CotasQHU(J,K1).. qhu(J,K1) =L= MHU(J)*YHUL(J,K1);

Equation CotasQCU(I,X) Cotas para la transferencia de calor de corrientes calientes con servicios de enfriamiento;
CotasQCU(I,K).. qcu(I,K) =L= MCU(I)*YCUL(I,K);


Equation VarBin(I,J);
VarBin(I,J).. Sum((K,K1), YL(I,K,J,K1)) =L= 1;

Equation VarBinHU(J);
VarBinHU(J).. Sum(K1, YHUL(J,K1)) =L= 1;

Equation VarBinCU(I,J);
VarBinCU(I,J).. Sum(K, YCUL(I,K)) =L= 1;


Equation TOTQHU; 
TOTQHU..  Sum((J,K1), qhu(J,K1)) =E= Sum((J,K1), qhu.l(J,K1)) ;


*Seccion de modelos
Model transportLP /objT, SUPPLY, DEMAND/;
Model transport /obj, SUPPLY, DEMAND, CotasQ, CotasQHU, CotasQCU/;
*Model transport /obj, SUPPLY, DEMAND, CotasQ, CotasQHU, CotasQCU, TOTQHU/;

*Solve transportLP using LP minimizing funobjL;
Solve transport using LP minimizing funobj;

