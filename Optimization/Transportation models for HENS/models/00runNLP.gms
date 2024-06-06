*-------------------------------------------------------------------------------
*--------------   GAMS FILE FOR TRANSPORTATION PROBLEM

*--------------   DEVELOPED BY LUIS E. PEDROZA ROBLES A.

*--------------   UNIVERSIDAD AUTÓNOMA METROPOLITANA IZTAPALAPA, MEXICO.
scalar starttime; starttime = jnow;
Scalar cumcputime;
Scalar flag, aux;


$inlinecom /* */
$eolcom <
$include options.c

Option NLP = CONOPT;

File PRUEBA /p-prueba.c/;
Scalar ATOT;
*----------------------------------------------------------------------------------
Scalar MimimumUtilities;   < 0 significa que no se consideran metas de energia
MimimumUtilities = 0;      < 1 significa que se resuelve el problema de uso mínimo de servicios y se fijan los valores de servicios


$include 01input11H-2C.dat


Parameter HRAT, DigHRAT, NO_ITERA;
option seed = 911;
HRAT =       10.000;
DigHRAT =  3;
NO_ITERA =   10;



*Probabilidad de que un apareamiento sea clonado sobre una corriente madre (0 indica que todos seran clonados, 1 indica que ninguno puede ser clonado)
Scalar ProBC /0/;
*----------------------------------------------------------------------------------

$include 02GenIntervalos.dat


*-------------------------------------------------------------------------------
*--------------      COEFICIENTES DE COSTO PARA EL PROBLEMA DE USO MÍNIMO DE SERTVICIOS AUXILIARES
Parameter CT(I,X,J,X1) Coeficientes de costo para recuperación de calor;
Parameter Chu(J,X1)   Coeficientes de costo para los servicios de calentamiento;
Parameter Ccu(I,X)    Coeficientes de costo para los servicios de enfriamiento;
Parameter BigM        /100000/;

Chu(J,K1) = 1;
Ccu(I,K) = 1;

CT(I,K,J,K1) = 0;
Loop[K,
     Loop[K1,
          If [OrdK(K) GT OrdK1(K1), CT(I,K,J,K1) = BigM];
             ];
     ];
*-------------

Scalar funobj;
Set SV /1*5000/;
set sol /multsol1*multsol100/; 
variables qhusol(sol,J,X1), qcusol(sol,I,X), qrecsol(sol,I,X,J,X1), funobjNLPsol(sol);

File TPFILE /P-transp.res/;
File METAS /P-Metas.gen/;

Parameter HTin Entalpía total que ingresa a la red
          HTout Entalpía total que sale de la red;

Parameter Res(X) Residuo para cada intervalo de temperatura
          Rfin Residuo final
          ITPP Intervalo de temperatura del punto de pliegue;

Parameter DTML(X,X1), DTMLH(X1), DTMLC(X);


Loop[(X,X1)$(Ord(X) LE (NOK) and Ord(X1) LE (NOK) and Ord(X1) GE Ord(X)),
  If((TH(X)-TC(X1)) EQ (TH(X+1)-TC(X1+1)),
    DTML(X,X1) = ((TH(X)-TC(X1))+(TH(X+1)-TC(X1+1)))/2;
  else 
    DTML(X,X1) = ((TH(X)-TC(X1))-(TH(X+1)-TC(X1+1)))/LOG((TH(X)-TC(X1))/(TH(X+1)-TC(X1+1)));
  );
];
DTML(K,K1)$(OrdK(K) GT OrdK1(K1)) = 1e-15;


Loop[I,
  Loop[X$(Ord(X) LE (NOK) and HS(I,X) GT 0),
    If(TH(X) GT TTCU and TH(X+1) GT TSCU,  
      If((TH(X)-TTCU) EQ (TH(X+1)-TSCU),
        DTMLC(X) = ((TH(X)-TTCU)+(TH(X+1)-TSCU))/2;
      Else 
        DTMLC(X) = ((TH(X)-TTCU)-(TH(X+1)-TSCU))/LOG((TH(X)-TTCU)/(TH(X+1)-TSCU));
      );
    Else
      If((TH(X)-TSCU) EQ (TH(X+1)-TSCU),
        DTMLC(X) = ((TH(X)-TSCU)+(TH(X+1)-TSCU))/2;
      Else
        DTMLC(X) = ((TH(X)-TSCU)-(TH(X+1)-TSCU))/LOG((TH(X)-TSCU)/(TH(X+1)-TSCU));
      );
    );
  ];
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


*-------------------------------------------------------------------------------
*--------------  MODELO DE TRANSPORTE
Parameter qhut Servicios totales de calentamiento
          qcut Servicios totales de enfriamiento
          qrec Calor total recuperado ;

Variable funobjNLP, MinUtil;


Positive variables
           q(I,X,J,X1) Calor recuperado
           qhu(J,X1)   Calor proveniente de servicios de calentamiento
           qcu(I,X)    Calor proveniente de servicios de enfriamiento;

Positive Variable A Area de recuperacion de calor
                  AHU Area de calentamiento
                  ACU Area de enfriamiento;


Positive variables NuHU, NuCU, NuRec;

Equation objU función a minimizar;
objU.. MinUtil =E= Sum((I,K,J,K1), CT(I,K,J,K1)*q(I,K,J,K1))
                +Sum((J,K1), Chu(J,K1)*qhu(J,K1))
                +Sum((I,K), Ccu(I,K)*qcu(I,K));


Equation objNLP función a minimizar;
objNLP..
  
funobjNLP =E= REC*(FANNUAL*CCFXCG+FANNUAL*CCPREC*A**CCAEXP)+
              NuHU*(FANNUAL*CCFXCGHU+FANNUAL*CCPRECHU*AHU**CCAEXPHU)+
              NuCU*(FANNUAL*CCFXCGCU+FANNUAL*CCPRECCU*ACU**CCAEXPCU)+
           Sum((J,K1)$(HD(J,K1) GT 0),HUCOST*qhu(J,K1))+Sum((I,K)$(HS(I,K) GT 0),CUCOST*qcu(I,K));  

Equation SUPPLY(I,X) Restricciones de suministro;
SUPPLY(I,K).. Sum((J,K1)$(OrdK1(K1) GE OrdK(K)),q(I,K,J,K1))+qcu(I,K) =E= HS(I,K);

Equation DEMAND(J,X1) Restricciones de demanda;
DEMAND(J,K1).. Sum((I,K)$(OrdK(K) LE OrdK1(K1)),q(I,K,J,K1))+qhu(J,K1) =E= HD(J,K1);

*Areas por equipo
Equation Area;
Area..  A =E= Sum((I,K,J,K1)$(HS(I,K) GT 0 and HD(J,K1) GT 0 and OrdK1(K1) GE OrdK(K)), 1/DTML(K,K1)*(1/FHTCI(I)+1/FHTCJ(J))*q(I,K,J,K1))/REC;

Equation AreaHU;
AreaHU.. AHU =E= Sum((J,K1)$(HD(J,K1) GT 0 and DTMLH(K1) GT 1e-5), 1/DTMLH(K1)*(1/FHTCHU+1/FHTCJ(J))*qhu(J,K1))/(NuHU+(1-NuHU/(NuHU+1e-6)));

Equation AreaCU;
AreaCU.. ACU =E= Sum((I,K)$(HS(I,K) GT 0 and DTMLC(K) GT 1e-5), 1/DTMLC(K)*(1/FHTCI(I)+1/FHTCCU)*qcu(I,K))/(NuCU+(1-NuCU/(NuCU+1e-6)));

Equation NEQHU;
NEQHU.. NuHU =E= Sum(J,Sum(K1, qhu(J,K1))/(Sum(K1,qhu(J,K1)+1e-6)));
*NEQHU.. NuHU =E= HEA;

Equation NEQCU;
NEQCU.. NuCU =E= Sum(I,Sum(K, qcu(I,K))/(Sum(K,qcu(I,K)+1e-6)));
*NEQCU.. NuCU =E= COO;


Equation qhumin;
qhumin.. sum((J,K1), qhu(J,K1)) =E= qhut;

Equation qcumin;
qcumin.. sum((I,K), qcu(I,K)) =E= qcut;


Model MinUtilLP /objU, SUPPLY, DEMAND/;
Model transportNLPUtil /ObjNLP, SUPPLY, DEMAND, Area, AreaHU, AreaCU, NEQHU, NEQCU, qhumin, qcumin/;
Model transportNLP /ObjNLP, SUPPLY, DEMAND, Area, AreaHU, AreaCU, NEQHU, NEQCU/;


q.lo(I,K,J,K1) = 0;
q.up(I,K,J,K1)$(OrdK1(K1) LT OrdK(K)) = 0;
q.up(I,K,J,K1)$(OrdK1(K1) GE OrdK(K)) = min(HS(I,K), HD(J,K1));

qhu.lo(J,K1) = 0;
qhu.up(J,K1) = HD(J,K1);

qcu.lo(I,K) = 0;
qcu.up(I,K) = HS(I,K);

NuHU.lo = 0.00000001;
NuHU.up = Card(J);

NuCU.lo = 0.00000001;
NuCU.up = Card(I);


Loop[SV$(Ord(SV) LE NO_ITERA),

  q.l(I,K,J,K1) = uniform(q.lo(I,K,J,K1), q.up(I,K,J,K1));
  qhu.l(J,K1) = uniform(qhu.lo(J,K1), qhu.up(J,K1));
  qcu.l(I,K) = uniform(qcu.lo(I,K), qcu.up(I,K));

  NuHU.l = uniform(NuHU.lo, NuHU.up);
  NuCU.l = uniform(NuCU.lo, NuCU.up);



*Areas por equipo
  A.l = Sum((I,K,J,K1), 1/DTML(K,K1)*(1/FHTCI(I)+1/FHTCJ(J))*q.l(I,K,J,K1))/REC;

  AHU.l = Sum((J,K1)$(DTMLH(K1) GT 1e-5), 1/DTMLH(K1)*(1/FHTCHU+1/FHTCJ(J))*qhu.l(J,K1))/(NuHU.l+(1-NuHU.l/(NuHU.l+1e-6))); 

  ACU.l = Sum((I,K)$(DTMLC(K) GT 1e-5), 1/DTMLC(K)*(1/FHTCI(I)+1/FHTCCU)*qcu.l(I,K))/(NuCU.l+(1-NuCU.l/(NuCU.l+1e-6)));



If (MimimumUtilities = 1,
  Solve MinUtilLP using LP minimizing MinUtil;


  qhut = Sum((J,K1), qhu.L(J,K1));
  qcut = Sum((I,K), qcu.L(I,K));
  qrec = Sum((I,K,J,K1), q.L(I,K,J,K1));


  transportNLPUtil.OptFile = 2;
  Solve transportNLPUtil using NLP minimizing funobjNLP;

Else 
  transportNLP.OptFile = 1;
  Solve transportNLP using NLP minimizing funobjNLP;

  );

funobj = funobjNLP.l;


qhut = Sum((J,K1), qhu.L(J,K1));
qcut = Sum((I,K), qcu.L(I,K));
qrec = Sum((I,K,J,K1), q.L(I,K,J,K1));



];

