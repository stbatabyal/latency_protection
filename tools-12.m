If[$VersionNumber < 6, <<Statistics`DataManipulation`];

(* list1 - list in which replacements with values in list2 needs to be done *)
(* pos indicates positions in which values in list1 will be replaced with *)
(* values in list2 *)

(* Does Arcsin/sqrt transformation for data given as percent *)
(* Inversing transformed data is done by InverseArcSinSqrt *)
(* Used in 
~/research/immunology/recirculation/lymph_nodes/Lymph_nodes-Tomura-2.nb *)

Clear[ImarisToCSV];

ColumnNames[tempdata_] := 
  Transpose[{Range[Length[tempdata]], tempdata}];

myHeader[dataset_] := 
  Transpose[{Range[Length[dataset[[1]]]], dataset[[1]]}];

ImarisToCSV[cellsfilename_, header_ : "Yes"] := 
  Module[{moviename, imfile, imfileData, imfileTime, 
    imfileAssociation, imfileTrackIDs, imfileTrackIDsLabels, labels, 
    times, tdata, tempdata1, tempdata2, posTrackID, posTrackIDtime}, 
   moviename = StringSplit[cellsfilename, "."][[1]];
   imfile = Import[cellsfilename, {"Sheets", "Position"}];
   imfileTime = Import[cellsfilename, {"Sheets", "Time"}];
   imfileAssociation = 
    AssociationThread[
     ToString /@ imfile[[2]] -> Range[Length[imfile[[2]]]]];
   imfileData = Rest[Rest[imfile]];
   imfileTrackIDs = 
    Union[imfileData[[All, imfileAssociation["TrackID"]]]];
   imfileTrackIDsLabels = 
    Map[moviename <> "." <> ToString[DecimalForm[Round[#]]] &, 
     imfileTrackIDs];
   
   posTrackID = Position[imfile[[2]], "TrackID"][[1, 1]];
   posTrackIDtime = Position[imfileTime[[2]], "TrackID"][[1, 1]];
   tdata = Table[
     tempdata1 = 
      Cases[imfile, x_?(#[[posTrackID]] == imfileTrackIDs[[track]] &)];
     tempdata2 = 
      Cases[imfileTime, 
       x_?(#[[posTrackIDtime]] == imfileTrackIDs[[track]] &)];
     times = tempdata1[[All, 7]];
     Table[
      Join[{imfileTrackIDsLabels[[track]]}, 
       Cases[tempdata1, 
         x_?(#[[7]] == times[[t]] &) :> Take[x, 3]][[1]], 
       Cases[tempdata2, x_?(#[[4]] == times[[t]] &) :> x[[1]]]], {t, 
       1, Length[times]}], {track, 1, Length[imfileTrackIDs]}];
   labels = {"Cell ID", "Time in " <> imfileTime[[3, 2]], 
     "Position X in " <> imfile[[3, 4]], 
     "Position Y in " <> imfile[[3, 4]], 
     "Position Z in " <> imfile[[3, 4]]};
   If[header == "Yes", Prepend[Flatten[tdata, 1], labels], 
    Flatten[tdata, 1]]
   ];


(* special case; κ = 0 is the uniform distribution *)
dimrothWatsonRandom[μ_?VectorQ, κ_ /; κ == 0] :=
                    Normalize[RandomVariate[NormalDistribution[], 3]]

(* https://mathematica.stackexchange.com/questions/13038/vectors-with-a-certain-magnitude-in-mathematica *)

(* splitting tracks for imaging data *)
(* dataset has 3-4 columns: 1st - is the trackID, then x,y,z,t or x,y,t *)
(* precision is the fraction by which time steps may differ between different time 
   slices. Typically, precision is 0.95 - i.e., time steps are nearly identical *)



Clear[cellSpeeds, cellTurningAngles, dataset, cellMSD, 
  cellVelocityCorrelation];


(* cell acceleration for every displacement *)

cellsAccelerations[dataset_, nd_ : 3] := 
  Map[If[Length[#] > 2, 
     Map[Sqrt[
         Total[(Take[#[[3]], nd] + Take[#[[1]], nd] - 
             2*Take[#[[2]], nd])^2]]/(#[[2, nd + 1]] - #[[1, 
             nd + 1]])^2 &, Partition[#, 3, 1]], "NA"] &, dataset];

(* cell speed per every displacement *)

cellSpeeds[dataset_, nd_ : 3] := 
  Map[If[Length[#] > 1, 
     Flatten[Map[
       EuclideanDistance[Take[#[[2]], nd], 
          Take[#[[1]], nd]]/(Take[#[[2]], -1] - Take[#[[1]], -1]) &, 
       Partition[#, 2, 1]]], {"NA"}] &, dataset];

(* turning angles per every displacement *)

cellTurningAngles[dataset_, nd_ : 3] := 
  Map[If[Length[#] > 2, 
     Map[180*VectorAngle[Take[#[[2]], nd] - Take[#[[1]], nd], 
          Take[#[[3]], nd] - Take[#[[2]], nd]]/Pi &, 
      Partition[#, 3, 1]], {"NA"}] &, dataset];


(* MSD for cohort *)

averageMSD[dataset_, nd_ : 3] := 
  Module[{alltimes, tmax, allpos, longesttrack, dt, tempdata2, 
    tempdata1, t1},
   alltimes = Map[Length[#] &, dataset];
   tmax = Max[alltimes];
   allpos = Position[alltimes, tmax];
   longesttrack = dataset[[allpos[[1, 1]]]];
   dt = Mean[
     Map[#[[2]] - #[[1]] &, 
      Partition[longesttrack[[All, nd + 1]], 2, 1]]];
   
   tempdata2 = Table[
     tempdata1 = dataset[[cell]];
     t1 = Length[tempdata1];
     Join[
      Table[{tempdata1[[t, -1]] - tempdata1[[1, -1]], 
        EuclideanDistance[Take[tempdata1[[1]], nd], 
          Take[tempdata1[[t]], nd]]^2}, {t, 2, t1}],
      Table[{i*dt, "NA"}, {i, t1 + 1, tmax}]
      ], {cell, 1, Length[dataset]}];
   Map[Mean[Cases[#, x_?(VectorQ[#, NumericQ] &)]] &, 
    Transpose[tempdata2]]
   ];



(* Calculating MSD per track *)

cellMSD[dataset_, nd_ : 3] := Module[{tempdata1,tempdata0},
  Table[
   tempdata0 = dataset[[cell]];
   Map[Mean[#] &,
    Split[Sort[
      Flatten[
       Table[
        tempdata1 = Drop[tempdata0, tt];
        Table[{tempdata1[[t, -1]] - tempdata1[[1, -1]], 
          EuclideanDistance[Take[tempdata1[[1]], nd], 
            Take[tempdata1[[t]], nd]]^2}, {t, 2, Length[tempdata1]}],
        {tt, 0, Length[tempdata0] - 2}], 1]], 
     Round[First[#1]] === Round[First[#2]] &]], {cell, 1, Length[dataset]}]
];

(* Calculating velocity correlation per track *)

cellVelocityCorrelation[dataset_, nd_ : 3] := Table[
   tempdata0 = dataset[[cell]];
   Prepend[Map[Mean[#] &,
     Split[Sort[Flatten[Table[
         tempdata1 = Drop[tempdata0, tt];
         tempdata2 = 
          Map[Append[
             Take[#[[2]], nd] - Take[#[[1]], nd], #[[2, -1]]] &, 
           Partition[tempdata1, 2, 1]];
         Table[{tempdata2[[t, -1]] - tempdata2[[1, -1]], 
           Cos[VectorAngle[Take[tempdata2[[1]], nd], 
             Take[tempdata2[[t]], nd]]]}, {t, 2, 
           Length[tempdata2]}], {tt, 0, Length[tempdata0] - 3}], 1]], 
      First[#1] === First[#2] &]
     ], {0, 1}], {cell, 1, Length[dataset]}
   ];


(* splitting tracks for imaging data *)
(* dataset has 3-4 columns: 1st - is the trackID, then x,y,z,t or x,y,t *)
(* precision is the fraction by which time steps may differ between different time 
   slices. Typically, precision is 0.95 - i.e., time steps are nearly identical *)

(* DT is the time step in the data, output - gives all timesteps per track or summary *)

MissingFramesTest[data0_, ndim_ : 3, DT_ : 1, output_ : 1] := 
 Module[{cells, datatemp1, datatemp2, datatemp3},
  cells = Split[Drop[data0[[All, 1]], 1]][[All, 1]];
  datatemp2 = Table[
    datatemp1 = 
     Cases[data0, x_?(#[[1]] == cells[[cell]] &) :> Drop[x, 1]];
    Length[
     Split[Sort[
        Map[Round[#[[2]] - #[[1]], DT] &, 
         Partition[datatemp1[[All, ndim + 1]], 2, 1]]]][[All, 1]]],
    {cell, 1, Length[cells]}
    ];
  Which[
   output == 1, 
   Transpose[{Range[1, Max[datatemp2]], 
     BinCounts[datatemp2, {1, Max[datatemp2] + 1}]}],
   output == 2, datatemp2
   ]
  ]

Clear[SplittingTracks];
(* data must be in format: trackID, x, y, (z), t *)
(* precision is how much timestep could vary *)
SplittingTracks[dataset_, precision_ : 0.95] := 
  Module[{dim = Length[dataset[[1]]] - 2, cells, datatemp0, datatemp1,
     datatemp2, tdiff, t, nc, newcells},
   datatemp0 = Split[dataset, #1[[1]] === #2[[1]] &];
   cells = Map[#[[1, 1]] &, datatemp0];
   datatemp2 =
    Table[
     datatemp1 = Map[Drop[#, 1] &, datatemp0[[cell]]];
     tdiff = 
      Median[Split[
         Sort[Map[#[[2]] - #[[1]] &, 
           Partition[datatemp1[[All, dim + 1]], 2, 1]]]][[All, 1]]];
     Clear[newcells, t, nc];
     t = 1;
     nc = 1;
     newcells[nc] = {datatemp1[[t]]};
     For[t = 2, t < Length[datatemp1], t++, 
      If[Abs[datatemp1[[t, dim + 1]] - datatemp1[[t - 1, dim + 1]]] >=
          tdiff*(1 - precision) && 
        Abs[datatemp1[[t, dim + 1]] - datatemp1[[t - 1, dim + 1]]] <= 
         tdiff*(1 + precision), 
       newcells[nc] = Append[newcells[nc], datatemp1[[t]]], nc++; (* 
       split the track *)
       newcells[nc] = {datatemp1[[t]]}]
      ];
     Flatten[
      Table[Map[
        Join[{ToString[DecimalForm[cells[[cell]]]] <> "." <> 
            ToString[nn]}, #] &, newcells[nn]], {nn, 1, nc}], 1],
     {cell, 1, Length[cells]}
     ];
   Flatten[datatemp2, 1]
   ];

(* this routine should be used if timestep in the movie is well-defined *)
SplittingTracks0[dataset_] :=
  Module[{dim = Length[dataset[[1]]] - 2, cells, datatemp0, datatemp1,
     dt, t, nc, newcells, allcleaneddataALL},
   dt = Min[
     Cases[Map[#[[2]] - #[[1]] &, 
       Partition[dataset[[All, -1]], 2, 1]], x_?(# > 0 &)]];
   datatemp0 = Split[dataset, #1[[1]] === #2[[1]] &];
   cells = Map[#[[1, 1]] &, datatemp0];
   timepos = dim + 1;
   allcleaneddataALL =
    Table[
     datatemp1 = Map[Drop[#, 1] &, datatemp0[[cell]]];
     Clear[newcells, t, nc];
     t = 1;
     nc = 1;
     newcells[nc] = {datatemp1[[t]]};
     For[t = 2, t < Length[datatemp1], t++,
      If[datatemp1[[t, timepos]] - datatemp1[[t - 1, timepos]] == dt,
       newcells[nc] = Append[newcells[nc], datatemp1[[t]]],
       nc++;
       newcells[nc] = {datatemp1[[t]]}
       ]
      ];
     Flatten[Table[
       Map[
        Join[{cells[[cell]] <> "." <> ToString[cell] <> "." <> 
            ToString[nn]}, #] &, newcells[nn]], {nn, 1, nc}], 1],
     {cell, 1, Length[cells]}
     ];
   Flatten[allcleaneddataALL, 1]
   ];


(* Timestep for each cell/track *)

AllTimeSteps[data0_] := Module[{cells, allDT},
   cells = Split[Drop[data0[[All, 1]], 1]][[All, 1]];
   allDT = 
    Table[Min[
      Map[#[[2]] - #[[1]] &, 
       Partition[
        Cases[data0, x_?(#[[1]] == cells[[cell]] &) :> x[[-1]]], 2, 
        1]]], {cell, 1, Length[cells]}];
   Transpose[{cells, allDT}]
   ];

vonMises2DRandom[\[Mu]_?VectorQ, \[Kappa]_?NumericQ] := 
  Block[{n\[Mu] = Normalize[\[Mu]], new\[Theta]},
   If[\[Kappa] >= 0, 
    new\[Theta] = RandomVariate[VonMisesDistribution[0, \[Kappa]]], 
    new\[Theta] = RandomVariate[VonMisesDistribution[0, -\[Kappa]]]];
   If[\[Kappa] >= 0, RotationTransform[new\[Theta], {0, 0}][n\[Mu]], 
    RotationTransform[new\[Theta], {0, 0}][-n\[Mu]]]
   ];

vonMisesFisherRandom[μ_?VectorQ, κ_?NumericQ] := Module[{ξ = RandomReal[], w},
  w = 1 + (Log[ξ] + Log[1 + (1 - ξ) Exp[-2 κ]/ξ])/κ;
  RotationTransform[{{0, 0, 1}, Normalize[μ]}][
   Append[Sqrt[1 - w^2] Normalize[RandomVariate[NormalDistribution[], 2]], w]]]


dimrothWatsonRandom[μ_?VectorQ, κ_ /; NumericQ[κ] && Positive[κ]] := Module[{c, u, v, w},
  c = Exp[-κ/2] Csch[κ/2]/2;
  While[
        {u, v} = RandomReal[1, 2]; w = Log[1 + u/c]/κ;
        v > Exp[κ w (w - 1)]];
  RotationTransform[{{0, 0, 1}, Normalize[μ]}][Append[
  Sqrt[1 - w^2] Normalize[RandomVariate[NormalDistribution[], 2]], RandomChoice[{-1, 1}] w]]]

dimrothWatsonRandom[μ_?VectorQ, κ_ /; NumericQ[κ] && Negative[κ]] := Module[{c, d, u, v, w},
  c = Sqrt[-κ]; d = ArcTan[c];
  While[
        {u, v} = RandomReal[1, 2]; w = Tan[d u]/c;
        v > (1 - κ w^2) Exp[κ w^2]];
  RotationTransform[{{0, 0, 1}, Normalize[μ]}][Append[
  Sqrt[1 - w^2] Normalize[RandomVariate[NormalDistribution[], 2]], RandomChoice[{-1, 1}] w]]]





MyMemberQ[list_, i_] := 
  Module[{ll = Cases[list, x_?(# == i &)]}, 
   If[Length[ll] > 0, True, False]];

(* Cluster size *)

clusters[dd_] := Table[{i, Length[Cases[dd, x_?(# == i &)]]}, {i, Min[dd], 
    Max[dd]}];
clusters0[dd_, {imin_, imax_}] :=  Table[{i, Length[Cases[dd, x_?(# == i &)]]}, {i, imin, imax}];
clustersN[dd_] := Table[{i, Length[Cases[dd, x_?(# == i &)]]/Length[dd]}, {i, Min[dd], Max[dd]}];



ArcSinSqrtOne[d_] := Map[{#[[1]], ArcSin[Sqrt[#[[2]]/100]]} &, d];
ArcSinSqrt[d_] := Map[ArcSinSqrtOne[#] &, d];

InverseArcSinSqrtOne[d_] := Map[{#[[1]], 100*(Sin[#[[2]]])^2} &, d];
InverseArcSinSqrt[d_] := Map[InverseArcSinSqrtOne[#] &, d];



(* Calculates CIs from bootstrap fits obtained using minimum2,3 functions *)

CalculateCIs[pars_, bestfit_, fits_, cis_] := Module[{tt, nb, i, bp,goodfits},
      goodfits = Cases[fits,x_?(NumericQ[#[[1]]]&)];
      tt = Map[pars /. #[[2]] &, goodfits];
      bp = pars /. bestfit[[2]];
      nb = Length[goodfits];
      Table[{pars[[i]],bp[[i]], 
          Sort[tt[[All, i]]][[Max[1, Round[nb*(1 - cis)/2]]]],
          Sort[tt[[All, i]]][[Max[1, Round[nb*(1 + cis)/2]]]]}, {i, 1,
          Length[pars]}]
      ];


(* Replaces values in the list1 by values in the list 2 in position pos *)

MyReplace[list1_, list2_, pos_] := Module[{i, init, final},
      If[Length[list2] == Length[pos],
        final = init = list1;
        For[i = 1, i <= Length[pos], i++, 
          final = ReplacePart[final, list2[[i]], pos[[i]]]],
        Print["List2 and pos are of different length"]
        ];
      final
      ];

Pvalue[sol_]:=sol[[1, 2, 1, 2, 4]];

ToLogOne[temp_] := Cases[temp,x_?(#[[2]]>0&):>{x[[1]], Log[10, x[[2]]]}];
ToLogOneX[temp_] := Cases[temp,x_?(#[[1]]>0&):>{Log[10,x[[1]]], x[[2]]}];
ToLogOneXY[temp_] := Cases[temp,x_?(#[[1]]>0&& #[[2]]>0 &):>{ Log[10,x[[1]]], Log[10, x[[2]]]}];

ToLog[temp_] := Map[ToLogOne[#] &, temp];
ToLogX[temp_] := Map[ToLogOneX[#] &, temp];
ToLogXY[temp_] := Map[ToLogOneXY[#] &, temp];

ImportCFSE[file_, nn_, par_] := Module[{data, n0, t0, tt1, tt2},
      data = Import[file, "Table"];
      data = If[par == 0, data, Transpose[data]];
      dataset = Delete[Transpose[Delete[data, 1]], 1];
      timedata = Delete[data[[1]], 1];
      n0 = Dimensions[dataset][[2]];
      t0 = Dimensions[dataset][[1]];
      Flatten[
        Table[{timedata[[i]], j, 
            If[j == nmax, Sum[dataset[[i, j + k]], {k, 1, n0 - nmax}], 
              dataset[[i, j + 1]]]}, {i, 1, t0}, {j, 0, nmax}], 1]
      ];

MyNumericQ[dataset_]:=Cases[dataset, _?(NumericQ[#] &)];

MyCleanOne[dataset_] := Cases[dataset, x_?(VectorQ[#, NumericQ] &)];

MyClean[dataset_] := Map[MyCleanOne[#]&,dataset];

MyDrop[list_,lab_,lab2_]:=Module[{aa},
aa = lab;
Table[aa = DeleteCases[aa, _?(# == lab2[[i]] &)], {i, 1, 
      Length[lab2]}];
MySelect[list,lab,aa]
];


MySelect0[list_,lab_,lab2_]:=Module[{ll},
ll = Flatten[Map[Position[lab, #] &, lab2]];
list[[ll]]
];

MySelect[list_, lab_, i_] := Module[{lab0},
    lab0=Map[ToString[#]&,lab];
    If[VectorQ[i],
      Table[
        If[MatrixQ[Position[lab0, ToString[i[[j]]]]], 
          First[Take[list, First[Position[lab0, ToString[i[[j]]]]]]], 
          "Null"], {j, 1, Length[i]}],
      If[MatrixQ[Position[lab0, ToString[i]]], 
        First[Take[list, First[Position[lab0, ToString[i]]]]], "Null"]
      ]
     ];

MySelect2[data_, lab_, sel_] := Module[{allpos},
    allpos = 
      Table[Position[Map[ToString[#] &, lab], ToString[sel[[i]]]][[1, 1]], {i,
           1, Length[sel]}];
    data[[allpos]]
    ];


(* Generates the boundaries if every subset has a particular length *)
permut[list_] := Table[Total[Take[list, i]], {i, 1, Length[list]}];


(* Generates list of boundaries (min,max) for each subclass *)
classes[list_] := 
    Table[{1, 0}*UnitStep[i - 1] + 
        Take[Drop[Prepend[permut[list], 1], i], 2], {i, 0, 
        Length[permut[list]] - 1}];


(* Subselects data for ranges given in the list range *)
(* Range can be generated using routine classes *)
(* range has to be in the form {{a1,a2},{a2,a3}...{an-1,an}} *)

AverageBinnedData::usage =
"Function provides the average y-values in the dataset assigning nbin (x,y) values per one bin. For example, nbin=10 will assign 10 values per each bin. The output is the average dataset with average x value per bin, average y value per bin and number of data points averaged";

AverageBinnedData[fitdata_, nbin_] := Module[{ranges, sdata},
      ranges = 
        Append[Table[nbin, {i, 1, Floor[Length[fitdata]/nbin]}], 
          Length[fitdata] - nbin*Floor[Length[fitdata]/nbin]];
      sdata = PartitionData[Sort[fitdata], ranges];
        Table[{Mean[sdata[[i]][[All, 1]]], 
            Total[sdata[[i]][[All, 2]]]/ranges[[i]], 
            Length[sdata[[i]][[All, 2]]]}, {i, 1, Length[ranges]}]
      ];

(* This adds estimate of the SD for the X variable *)
AverageBinnedDataXY[fitdata_, nbin_] := Module[{ranges, sdata},
      ranges =
        Append[Table[nbin, {i, 1, Floor[Length[fitdata]/nbin]}],
          Length[fitdata] - nbin*Floor[Length[fitdata]/nbin]];
      sdata = PartitionData[Sort[fitdata], ranges];
        Table[{
            {Mean[sdata[[i]][[All, 1]]],
     	     StandardDeviation[sdata[[i]][[All, 1]]]},
            {Total[sdata[[i]][[All, 2]]]/ranges[[i]],
            Length[sdata[[i]][[All, 2]]]}
            }, {i, 1, Length[ranges]}]
      ];

AverageBinnedDataEqual::usage =
"Function provides the average y-values in the dataset when all x-values are divided in nbin intervals of the same length. It uses function SelectDataN.";

AverageBinnedDataEqual[fitdata_, nbin_] := Module[{tt, xmin, xmax},
      tt = fitdata[[All, 1]];
      xmin = Floor[Min[tt]];
      xmax = Ceiling[Max[tt]];
      tt=SelectDataN[fitdata, 
        Table[{i, i + (xmax - xmin)/nbin}, {i, xmin, 
            xmax - (xmax - xmin)/nbin, (xmax - xmin)/nbin}]];
      Cases[tt,_?(VectorQ[#,NumericQ]&)]
      ];

SelectDataN::usage = "Function provides the average y-values in the dataset assigning ranges for data to be binned. Ranges have to be specified as set of lists {{x1,x2},{x2,x3},{x3,x4}...}. The data are averaged per given x-interval";

SelectData[data_, range_] := 
    Table[{Mean[range[[i]]], 
        Mean[Cases[data, 
              x_?(#[[1]] >= range[[i, 1]] && #[[1]] <= range[[i, 2]] &) :> 
                x][[All, 2]]]}, {i, 1, Length[range]}];

(* This adds the number of elements per category *)
SelectDataN[data_, range_] := Module[{tt},
    Table[{Mean[range[[i]]], 
        tt=Cases[data, 
              x_?(#[[1]] >= range[[i, 1]] && #[[1]] <= range[[i, 2]] &) :> 
                x][[All, 2]];
         Mean[tt],Length[tt]}, {i, 1, Length[range]}]
     ];


(* This has an estimate for the average X and +/- SD *)

SelectDataNXY[data_, range_] := Module[{tt},
    Table[
	tt=Cases[data,
              x_?(#[[1]] >= range[[i, 1]] && #[[1]] <= range[[i, 2]] &)];
        {
         {Mean[tt[[All,1]]],StandardDeviation[tt[[All,1]]]},
         {Mean[tt[[All,2]]],Length[tt]}
        }, 
       {i, 1, Length[range]}]
     ];

(* data have to be in the form (number.mut, number.total.sequences *)

CI=0.95;
range[temp_] := Module[{q1, q2, x, n},
      {x, n} = temp;
      q2 = {If[x == 0, 0, 
            If[x == n, ((1 - CI)/2)^(1/n), 
              Quantile[BetaDistribution[x + 1/2, n - x + 1/2], (1 - CI)/2]]], 
          If[x == 0, 1 - ((1 - CI)/2)^(1/n), 
            If[x == n, 1, 
              Quantile[
                BetaDistribution[x + 1/2, n - x + 1/2], (1 + CI)/2]]]};
      q1 = n*q2;
      {q1, q2}
      ];

(* Add binomial error bars to the data *)

AddErrorBar[datatemp_] := 
  Map[{#[[1]], 
     Around[#[[2]], Last[range[{#[[2]]*#[[3]], #[[3]]}]] - #[[2]]]} &,
    datatemp];

AddErrorBarLog[datatemp_] :=
  Map[{#[[1]],
     Around[Log10[#[[2]]], Log10[Last[range[{#[[2]]*#[[3]], #[[3]]}]]] - Log10[#[[2]]]]} &,
    datatemp];


(* Add binomial error bars to the Y data and Normal error bars to X data*)

AddErrorBarXY[datatemp_] := 
  Map[{Around[#[[All, 1]][[1]], 
     {-If[CI == 0.95, 2, 1]*#[[1, 2]], 
       If[CI == 0.95, 2, 1]*#[[1, 2]]}],
       Around[#[[All, 1]][[2]], 
      Last[range[{#[[2, 1]]*#[[2, 2]], #[[2, 2]]}]] - #[[2, 1]]]} &, 
   datatemp];

(* Export data in table format to print in prism or similar *) 

AddErrorBarXYPrism95[datatemp_] :=
 Prepend[(Flatten[{#1[[All, 1]][[1]], #1[[All, 1]][[1]] - 
        If[CI == 0.95`, 2, 1] #1[[1, 2]], #1[[All, 1]][[1]] + 
        If[CI == 0.95`, 2, 1] #1[[1, 2]], #1[[All, 1]][[2]], 
       Last[range[{#1[[2, 1]] #1[[2, 2]], #1[[2, 2]]}]]}] &) /@ 
   datatemp, {"x value", "x 95% low", "x 95% high", "y value", 
   "y 95% low", "y 95% high"}];

AddErrorBarPrism[datatemp_] :=
    Prepend[
     Map[Join[Drop[#, -1],
          Last[range[{#[[2]]*#[[3]], #[[3]]}]]] &,
      datatemp], {"x value", "y value", "y lower limit",
    "y upper limit"}];

AddErrorBarXYPrism[datatemp_] := 
  Prepend[(Flatten[{#1[[All, 1]][[1]], If[CI==0.95,2,1]*#1[[1, 2]], 
        #1[[All, 1]][[2]], 
        Last[range[{#1[[2, 1]] #1[[2, 2]], #1[[2, 2]]}]]}] &) /@ 
    datatemp, {"x value", "Err.bar", "y value", "y lower limit", 
    "y upper limit"}];


AverageBinnedDataXY0[fitdata_, nbin_] := 
 Module[{ranges, sdata}, 
  ranges = Append[Table[nbin, {i, 1, Floor[Length[fitdata]/nbin]}], 
    Length[fitdata] - nbin Floor[Length[fitdata]/nbin]]; 
  sdata = PartitionData[Sort[fitdata], ranges]; 
  Table[{{Mean[sdata[[i]][[All, 1]]], 
     StandardDeviation[sdata[[i]][[All, 1]]]}, {Mean[
      sdata[[i]][[All, 2]]], 
     StandardDeviation[sdata[[i]][[All, 2]]]}}, {i, 1, 
    Length[ranges]}]]



AddErrorBarXY0[datatemp_] := Map[{
     #[[All, 1]], 
     ErrorBar[{-If[CI == 0.95, 2, 1]*#[[1, 2]], 
       If[CI == 0.95, 2, 1]*#[[1, 2]]},
      {-If[CI == 0.95, 2, 1]*#[[2, 2]], If[CI == 0.95, 2, 1]*#[[2, 2]]}
      ]
     } &, datatemp];

AddErrorBarXY0[datatemp_] := Map[Map[Around[#[[1]], 
     If[CI == 0.95, 2, 1]*#[[2]]] &, #] &, datatemp]


Clear[column,columnDrop];
If[$VersionNumber<6,
   column[x_,y_]:=Column[x,y];
   columnDrop[x_,y_]:=ColumnDrop[x,y], 
   column[x_,y_]:=Map[#[[y]]&,x];
   columnDrop[x_,y_]:=x[[All,Delete[Range[Length[x[[1]]]],y]]]
];

MedianDataOne[list_] := Module[{times, ntimes},
      times = column[Split[Sort[column[list, 1]]], 1];
      ntimes = Length[times];
      Table[{times[[i]], 
          Median[Cases[list, x_?(#[[1]] == times[[i]] &) :> x[[2]]]]}, {i, 1, 
          ntimes}]
      ];

AverageDataOneN[list_] := 
  Module[{times, ntimes}, 
   times = column[Split[Sort[column[list, 1]]], 1];
   ntimes = Length[times];
   Table[{times[[i]], 
     Mean[Cases[list, x_?(#[[1]] == times[[i]] &) :> x[[2]]]], 
     Length[Cases[list, x_?(#[[1]] == times[[i]] &) :> x[[2]]]]}, {i, 
     1, ntimes}]];

AverageDataOne[list_] := Module[{times, ntimes},
      times = column[Split[Sort[column[list, 1]]], 1];
      ntimes = Length[times];
      Table[{times[[i]], 
          Mean[Cases[list, x_?(#[[1]] == times[[i]] &) :> x[[2]]]]}, {i, 1, 
          ntimes}]
      ];

SDDataOne[list_] := Module[{times, ntimes},
      times = column[Split[Sort[column[list, 1]]], 1];
      ntimes = Length[times];
      Table[{times[[i]],
          StandardDeviation[Cases[list, x_?(#[[1]] == times[[i]] &) :> x[[2]]]]}, 
          {i, 1, ntimes}]
      ];

MinMaxDataOne[list_] := Module[{times, ntimes},
      times = column[Split[Sort[column[list, 1]]], 1];
      ntimes = Length[times];
      Table[{
          times[[i]],
          Min[Cases[list, x_?(#[[1]] == times[[i]] &) :> x[[2]]]],
          Max[Cases[list, x_?(#[[1]] == times[[i]] &) :> x[[2]]]]
          },
          {i, 1, ntimes}]
      ];


AverageData[list_] := Map[AverageDataOne[#] &, list];

MedianData[list_] := Map[MedianDataOne[#] &, list];

SDData[list_]:=Module[{times,ntimes},
times = column[Split[Sort[column[list[[1]], 1]]], 1];
ntimes = Length[times];
Map[Table[{times[[i]], StandardDeviation[Cases[#, x_?(#[[1]] == times[[i]] &) :> x[[2]]]]},
{i, 1, ntimes}] &, list]
];



(* Partitions data into sublists of length specified in divs 
    
    divs = {5,5,5} will create 3 sublists each of length 5

*)
PartitionData[temp_,divs_]:=partdata[temp,divs];

partdata[temp_, divs_] := 
    Table[Take[temp, classes[divs][[i]]], {i, 1, Length[classes[divs]]}];


ToLatex[list_, norm_, coef___] := Module[{coef0},
      If[Length[coef] == 0, coef0 = Table[1, {Length[list]}], coef0 = 1/coef];
      If[Length[list] != Length[norm],
        Print["Not Equal Lenghts"],
        Table[{list[[i, 1]], 
            MyRound[coef0[[i]]*list[[i, 2]], norm[[i]]], "(" <> 
              ToString[MyRound[coef0[[i]]*list[[i, 3]], norm[[i]]]] <> "--" <>
               ToString[MyRound[coef0[[i]]*list[[i, 4]], norm[[i]]]]<>")"}, {i, 1, 
            Length[list]}]
        ]
      ];


ToLatex2[list_, norm_, coef___] := Module[{coef0},
      If[Length[coef] == 0, coef0 = Table[1, {Length[list]}], coef0 = 1/coef];
      If[Length[list] != Length[norm],
        Print["Not Equal Lenghts"],
        Table[{ToString[list[[i, 1]]]<>"*"<>ToString[coef[[i]]], 
            ToString[MyRound[coef0[[i]]*list[[i, 2]], norm[[i]]]] <> " (" <> 
              ToString[MyRound[coef0[[i]]*list[[i, 3]], norm[[i]]]] <> "--" <>
               ToString[MyRound[coef0[[i]]*list[[i, 4]], norm[[i]]]]<>")"}, 
                {i, 1, 
            Length[list]}]
        ]
      ];


ToLatexNoLabel[list_, norm_,coef___] := Module[{coef0},
    If[Length[coef] == 0, coef0 = Table[1, {Length[list]}], coef0 = 1/coef];
    If[Length[list] != Length[norm], Print["Not Equal Lenghts"]];
    Table[{ToString[MyRound[coef0[[i]]*list[[i, 1]], norm[[i]]]] <> " (" <>
              ToString[MyRound[coef0[[i]]*list[[i, 2]], norm[[i]]]] <> "--" <>
               ToString[MyRound[coef0[[i]]*list[[i, 3]], norm[[i]]]]<>")"},
                {i, 1, Length[list]}]
        ];

AddColumns[set_,il1_]:=If[!MatrixQ[il1],
                Transpose[Append[Transpose[set], il1]],
                Transpose[Flatten[{Transpose[set], Transpose[il1]}, 1]]];

MyRound[x_,y_]:=N[Round[10^y*x]/10^y]


(* Steady state of the system *)

SteadyState[eqns_,vars_]:= Solve[Table[eqns[[i, 2]] == 0, {i, 1, Length[eqns]}], vars];
