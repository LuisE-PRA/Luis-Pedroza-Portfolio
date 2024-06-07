*-------------------------------------------------------------------------------
*--------------   GAMS FILE FOR TABU SEARCH

*--------------   DEVELOPED BY LUIS E. PEDROZA ROBLES A.

*--------------   UNIVERSIDAD AUTÓNOMA METROPOLITANA IZTAPALAPA, MEXICO.
$include 01input20SP1.dat
$include 02GenIntervalos.dat

Option LP = Cplex;
*Option LP = Xpress;
*Option LP = Minos;

*-------------------------------------------------------------------------------
*--------------      COEFICIENTES DE COSTO PARA EL PROBLEMA DE TRANSPORTE
Parameter CT(I,X,J,X1) Coeficientes de costo para recuperación de calor;
Parameter Chu(J,X1)   Coeficientes de costo para los servicios de calentamiento;
Parameter Ccu(I,X)    Coeficientes de costo para los servicios de enfriamiento;
Parameter BigM        /100000/;

Chu(J,K1) = 1;
Ccu(I,K) = 1;

CT(I,K,J,K1) = 0;
Loop[K,
     Loop[K1,
          If [OrdK(K) GT OrdK(K1), CT(I,K,J,K1) = BigM];
             ];
     ];



*-------------------------------------------------------------------------------
*--------------  MODELO DE TRANSPORTE
Variable funobj;
Positive variables
           q(I,X,J,X1) Calor recuperado
           qhu(J,X1)   Calor proveniente de servicios de calentamiento
           qcu(I,X)    Calor proveniente de servicios de enfriamiento;

Equation obj función a minimizar;
obj..
funobj =E= Sum((I,K,J,K1),
                         CT(I,K,J,K1)*q(I,K,J,K1))
          +Sum((J,K1), Chu(J,K1)*qhu(J,K1))
          +Sum((I,K), Ccu(I,K)*qcu(I,K));


Equation SUPPLY(I,X) Restricciones de suministro;
SUPPLY(I,K).. Sum((J,K1),q(I,K,J,K1))+qcu(I,K) =E= HS(I,K);

Equation DEMAND(J,X1) Restricciones de demanda;
DEMAND(J,K1).. Sum((I,K),q(I,K,J,K1))+qhu(J,K1) =E= HD(J,K1);

$include options.c

Model transport /obj,SUPPLY,DEMAND/;

Solve transport using LP minimizing funobj;




*-------------------------------------------------------------------------------
*--------------      OUTPUT FILE FOR TRANSPORT PROBLEM

File TPFILE /z-transp.res/;
Put  TPFILE;


*-------------------------------------------------------------------------------
*---------        PROCESAMIENTO DE DATOS

Put ///
$If exist z-Datos.c Put @30 'HRAT = ', HRAT:12:8 //;
Put      #3,@27  'UNIVERSIDAD AUTÓNOMA METROPOLITANA'
         #5,@42  'Proyecto'
         #7,@35  'Caso de estudio 10SP1'
         #9,@21  'Autores: Luis E. Pedroza, Juan M. Zamora Mata'

         #14,@1'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-'/
         #15,@38'Archivo de resultados'
         #16,@1'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-'///;

Put @10 'Información de corrientes y costos para el caso de estudio'//

@1 'Corriente', @12 'TS (°C)', @22 'TO (°C)', @32 'FCp (kW/K)', @45 'h (kW/Km^2)', @ 59 'qdisp/qrec (kW)'/
'_________________________________________________________________________'/;

Loop[I,
    Put @5 'H', Ord(I):1:0, @14 TSI(I):3:0, @23 TTI(I):3:0, @33 FCPI(I):6:2, @47 FHTCI(I):5:3, @61 QTOTI(I):10:3
     /];
Loop[J,
    Put @5 'C', Ord(J):1:0, @14 TSJ(J):3:0, @23 TTJ(J):3:0, @33 FCPJ(J):6:2, @47 FHTCJ(J):5:3, @61 QTOTJ(J):10:3
     /];
Put/;
Put @5 'HU', @14 TSHU:3:0, @23 TTHU:3:0, @47 FHTCHU:5:3/
    @5 'CU', @14 TSCU:3:0, @23 TTCU:3:0, @47 FHTCCU:5:3/
'_________________________________________________________________________'///;

Parameter qhut Servicios totales de calentamiento
          qcut Servicios totales de enfriamiento
          qrec Calor total recuperado ;

qhut = Sum((J,K1), qhu.L(J,K1));
qcut = Sum((I,K), qcu.L(I,K));
qrec = Sum((I,K,J,K1), q.L(I,K,J,K1));

Put '----------------Metas de energía--------------------'//
    qrec.ts ' = ', qrec:15:3 /
    qhut.ts ' = ', qhut:15:3/
    qcut.ts ' = ', qcut:15:3///;
*-------------------------------------------------------------------------------
Parameter HTin Entalpía total que ingresa a la red
          HTout Entalpía total que sale de la red;

HTin = Sum[(I,K), HS(I,K)]+qhut;
HTout = Sum[(J,K1), HD(J,K1)]+qcut;

Put '------------Balance global de energía---------------'//
    HTin.ts ' = ', HTin:15:3/
    HTout.ts ' = ', HTout:15:3;
*-------------------------------------------------------------------------------
Parameter Res(X) Residuo para cada intervalo de temperatura
          Rfin Residuo final
          ITPP Intervalo de temperatura del punto de pliegue;


Res('1') = qhut+Sum[I, HS(I,'1')]-Sum[J, HD(J,'1')];

Loop[X$(Ord(X) GT 1 and (Ord(X) LE NOK)),
Res(X) = Res(X-1)+Sum[I, HS(I,X)]-Sum[J, HD(J,X)];
      ];

Put ///'----------------Punto de Pliegue--------------------'//;
Put @1'Intervalo', @20'Residuo'/;

Loop[K,
Put OrdK(K):5:0, Put @15 Res(K)/;
     ];

Loop[K,
      Rfin = Res(K)-qcut];

Put / Rfin.ts " = ",  Rfin:4:2;
Put //;
Put 'El punto de pliegue se ubica entre los intervalos de temperatura ';

Loop[K,
        If [Res(K) LE EPS,
            ITPP = OrdK(K)
            Put ITPP:1:0 ' y '  Put (ITPP+1):1:0 /;
            ];
     ];


Put 'En las Temperaturas [';
Loop [X$(Ord(X) LE NOK),
      If[Ord(X) EQ ITPP,
         Put TH(X+1):8:4,', ' TC(X+1):8:4;
         ];
      ];
Put ']'///;

Put '--------------Tabla 1. Apareaminetos sugeridos--------------------'//
'Corriente caliente I desde intervalo K intercambia'/
'calor con corriente fría J en intervalo K1'//

     @3 'I', @6 'K', @9 'J', @11 'K1', @22 'q'/;

Loop[J,
     Loop[K1,
          If[qhu.L(J,K1) <> 0,
                           Put @1 '-qhu-', @8 Ord(J):2:0,
                               @11 OrdK(K1):2:0, @16 qhu.L(J,K1):10:3/];
          ];
     ];
Put '-------------------------'/;

*Si no se desea distinguir el PP
Loop [I,
      Loop [J,
            Loop [K,
                  Loop [K1,
                        If[q.L(I,K,J,K1) <> 0,
                           Put @2 Ord(I):2:0,@5 OrdK(K):2:0, @8 Ord(J):2:0,
                               @11 OrdK(K1):2:0, @16 q.L(I,K,J,K1):10:3/];
                        ];
                  ];
            ];
      ];


$ontext
*Si se desea duistinguir el Punto de Pliegue
Loop [I,
      Loop [J,
            Loop [K$(OrdK(K) LE NOK),
                  Loop [K1$(OrdK(K1) LE NOK),
                        If[q.L(I,K,J,K1) <> 0 and OrdK(K) LE ITPP and OrdK(K1) LE ITPP,
                           Put @2 Ord(I):2:0,@5 OrdK(K):2:0, @8 Ord(J):2:0,
                               @11 OrdK(K1):2:0, @16 q.L(I,K,J,K1):10:3/];
                        ];
                  ];
            ];
      ];

Put '------------PP-----------'/;
Loop [I,
      Loop [J,
            Loop [K$(OrdK(K) LE NOK),
                  Loop [K1$(OrdK(K1) LE NOK),
                        If[q.L(I,K,J,K1) <> 0 and OrdK(K) GT ITPP and OrdK(K1) GT ITPP,
                           Put @2 Ord(I):2:0,@5 OrdK(K):2:0, @8 Ord(J):2:0,
                               @11 OrdK(K1):2:0, @16 q.L(I,K,J,K1):10:3/];
                        ];
                  ];
            ];
      ];

$Offtext
Put '-------------------------'/;
Loop[I,
     Loop[K,
          If[qcu.L(I,K) <> 0,
                           Put @2 Ord(I):2:0,@5 OrdK(K):2:0,
                           @8 '-qcu-', @16 qcu.L(I,K):10:3/];
          ];
     ];

*Si se requieren servicios extremos en todas las corrientes
If [OPSERV = 2, qhu.L(J,K1) = 1; qcu.L(I,K) = 1];

