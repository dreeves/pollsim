#!/usr/bin/env mash

Needs["MultivariateStatistics`"];
(* Needs["Histograms`"]; *)

pathfunc[_] = "./";
pathfunc["yootles"] = "/var/www/html/yootles/pollsim/";
pathfunc["danny"] = pathfunc["dreev"] = "/Users/dreeves/prj/pollsim/";
path = pathfunc[$MachineName];
xpath = "yootles.com:/var/www/yootles/pollsim";

EDAY = {2012,11,6};  (* election day *)
NSIM = 4000; (* number of simulations to run *)
PMIN = .3;  (* intrade price below which true prob is zero; originally .3 *)
PMAX = .7;  (* intrade price above which true prob is one; originally .7 *)
WIND = 90;
  
ptsz = .007;  (* point size for plotted points in graph *)
  
mndist = MultinormalDistribution;
  
(* Intrade contract IDs for democrat and republican victories by state *)
(*
states = 
  {"AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", 
   "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", 
   "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", 
   "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", 
   "WV", "WI", "WY"};
itd = 
  {745711, 745714, 745717, 745720, 745723, 745726, 745729, 766056, 745732,
   745735, 745738, 745741, 745744, 745747, 745750, 745753, 745756, 745759,
   745762, 745765, 745768, 745771, 745774, 745777, 745780, 745783, 745786,
   745789, 745792, 745795, 745798, 745801, 745804, 745807, 745810, 745813,
   745816, 745819, 745822, 745825, 745828, 745831, 745834, 745837, 745840,
   745843, 745846, 745849, 745852, 745855, 745858};
itr =
  {745712, 745715, 745718, 745721, 745724, 745727, 745730, 766057, 745733,
   745736, 745739, 745742, 745745, 745748, 745751, 745754, 745757, 745760,
   745763, 745766, 745769, 745772, 745775, 745778, 745781, 745784, 745787,
   745790, 745793, 745796, 745799, 745802, 745805, 745808, 745811, 745814,
   745817, 745820, 745823, 745826, 745829, 745832, 745835, 745838, 745841,
   745844, 745847, 745850, 745853, 745856, 745859}; 
*)
(* from 2008:
itd = 
  {416468, 416471, 416484, 416491, 416498, 416505, 416508, 416623, 416511, 
   417861, 416514, 416520, 416523, 416532, 416542, 416545, 416554, 416557, 
   416560, 416569, 416584, 416593, 416599, 416602, 416612, 416619, 416487, 
   416490, 416496, 416502, 416516, 416526, 416529, 416534, 416539, 416548, 
   416551, 416590, 416595, 416605, 416608, 416611, 416617, 416632, 416636, 
   416639, 416642, 416645, 416648, 416653, 416654};
itr =
  {416469, 416472, 416485, 416493, 416500, 416506, 416509, 416624, 416512, 
   417866, 416515, 416521, 416524, 416537, 416543, 416546, 416555, 416558, 
   416561, 416570, 416585, 416594, 416600, 416603, 416614, 416621, 416488, 
   416492, 416497, 416503, 416518, 416527, 416530, 416535, 416540, 416549, 
   416552, 416591, 416597, 416606, 416609, 416613, 416618, 416634, 416637, 
   416640, 416643, 416646, 416649, 416651, 416655};
*)

(* Electoral College events: *)
RepECID = Range[756126,756146]; (* 180 to 380; shouldn't hardcode that! *)
(* from 2008:
RepECID = {648314, 648313, 648312, 613052, 613053, 613054, 613055, 613056, 
           613057, 613058, 613059, 613060, 613061, 613062, 613777, 613778, 
           613779, 613780, 613781, 613782, 613783};
*)
DemECID = Range[756101,756120]; (* 210 to 400 electoral votes *)
(* from 2008:
DemECID = {613041, 613042, 613043, 613044, 613045, 613046, 613047, 613048, 
           613049, 613050, 613051, 613770, 613771, 613772, 613773, 613774, 
           613775, 613776, 649913, 649914};
*)
   
(* Number of electoral votes by state. Fetched from here: *)
(* evUrl = "http://www.electoral-vote.com/evp2008/Pres/Excel/today.csv"; *)
(*
ev = {9, 3, 10, 6, 55, 9, 7, 3, 3, 27, 15, 4, 4, 21, 11, 7, 6, 8, 9, 4, 10, 
      12, 17, 10, 6, 11, 3, 5, 5, 4, 15, 5, 31, 15, 3, 20, 7, 7, 21, 4, 8, 
      3, 11, 34, 5, 3, 13, 11, 5, 10, 3};
*)

(* State, Intrade ID for dem winner, Ditto for R, Electoral votes, Shoo-in *)
(* pasted from a spreadsheet from cody emailed on 2012-10-08 *)
data2012 = ImportString[
  "AL\t745711\t745712\t9\tR
   AK\t745714\t745715\t3\tR
   AZ\t745717\t745718\t11\t
   AR\t745720\t745721\t6\tR
   CA\t745723\t745724\t55\tD
   CO\t745726\t745727\t9\t
   CT\t745729\t745730\t7\tD
   DC\t766056\t766057\t3\tD
   DE\t745732\t745733\t3\tD
   FL\t745735\t745736\t29\t
   GA\t745738\t745739\t16\tR
   HI\t745741\t745742\t4\tD
   ID\t745744\t745745\t4\tR
   IL\t745747\t745748\t20\tD
   IN\t745750\t745751\t11\t
   IA\t745753\t745754\t6\t
   KS\t745756\t745757\t6\tR
   KY\t745759\t745760\t8\tR
   LA\t745762\t745763\t8\tR
   ME\t745765\t745766\t4\tD
   MD\t745768\t745769\t10\tD
   MA\t745771\t745772\t11\tD
   MI\t745774\t745775\t16\tD
   MN\t745777\t745778\t10\tD
   MS\t745780\t745781\t6\tR
   MO\t745783\t745784\t10\t
   MT\t745786\t745787\t3\t
   NE\t745789\t745790\t5\tR
   NV\t745792\t745793\t6\t
   NH\t745795\t745796\t4\t
   NJ\t745798\t745799\t14\tD
   NM\t745801\t745802\t5\tD
   NY\t745804\t745805\t29\tD
   NC\t745807\t745808\t15\t
   ND\t745810\t745811\t3\tR
   OH\t745813\t745814\t18\t
   OK\t745816\t745817\t7\tR
   OR\t745819\t745820\t7\tD
   PA\t745822\t745823\t20\tD
   RI\t745825\t745826\t4\tD
   SC\t745828\t745829\t9\tR
   SD\t745831\t745832\t3\tR
   TN\t745834\t745835\t11\tR
   TX\t745837\t745838\t38\tR
   UT\t745840\t745841\t6\tR
   VT\t745843\t745844\t3\tD
   VA\t745846\t745847\t13\t
   WA\t745849\t745850\t12\tD
   WV\t745852\t745853\t5\tR
   WI\t745855\t745856\t10\t
   WY\t745858\t745859\t3\tR", "Table"];

states = data2012[[All,1]];
itd    = data2012[[All,2]];
itr    = data2012[[All,3]];
ev     = data2012[[All,4]];

(* Parse the above so that shoo[state] gives "D" or "R" or Null depending on if
   the democrat or republican or neither is a shoo-in in that state.
   Similarly, shoop[id] gives what the price would be for a contract id if there
   were any trading (there's no trading because D or R is a shoo-in).
   Also make a hash from ID back to state abbrevation. *)
each[x_, data2012, 
  statehash[x[[2]]] = x[[1]];
  statehash[x[[3]]] = x[[1]];
  partyhash[x[[2]]] = "D";
  partyhash[x[[3]]] = "R";
  Which[
    Length[x]<=4, shoo[x[[1]]] = Null,
    x[[5]]=="D",  shoo[x[[1]]] = "D"; shoop[x[[2]]] = 100; shoop[x[[3]]] = 0,
    x[[5]]=="R",  shoo[x[[1]]] = "R"; shoop[x[[2]]] = 0;   shoop[x[[3]]] = 100]]

perturb[p_] := Which[
  (* True,   p, *)
  p==100, RandomReal[{99.995,99.999}],
  p==0,   RandomReal[{0.0005,0.001}],
  True,   RandomReal[{p-0.0001,p+0.0001}]
]

fakedata[p_] := With[{start = {2012,2,7}},
  Table[{DateString[DatePlus[start, i], 
                    {"MonthNameShort", " ", "DayShort", ", ", "Year"}],
         perturb@p, perturb@p, perturb@p, perturb@p, 666}, 
        {i, 0, DateDifference[start, DateList[]]}]]

(*************************************************************************)
(*                               INTRADE API                             *)
(*************************************************************************)
    
(* intrade (price) history URL *)
ithUrl = 
 "http://data.intrade.com/graphing/jsp/downloadClosingPrice.jsp?contractId=";

(* Get historical prices. Rows are {Date, Open, Low, High, Close, Volume}. *)
phist[id_] := Module[{p = partyhash[id], s = shoo[statehash[id]]},
  Which[
    p == "D" && s == "D", Return@fakedata[100],
    p == "D" && s == "R", Return@fakedata[0],
    p == "R" && s == "D", Return@fakedata[0],
    p == "R" && s == "R", Return@fakedata[100]];
  Check[Rest@Import[cat[ithUrl, id]], 
    prn["\n\n<h1>phist: no intrade data; check back in 10 minutes</h1>\n"];
    Exit[1]]]

(* Get electoral college threshold market prices. *)
echist[id_] := Check[Import[cat[ithUrl, id]],
  prn["\n\n<h1>echist: no intrade data; check back in 10 minutes</h1>\n"];
  Exit[1]]

(* intrade contract info URL. cf: intrade.com/aav2/api/MarketMaker.pdf *)
itcUrl = "http://www.intrade.com/jsp/XML/MarketData/ConInfo.jsp?id=";
(* intrade price info URL *)
itpUrl = cat["http://www.intrade.com/jsp/XML/MarketData/ContractBookXML.jsp?",
             "depth=1&timestamp=0&id="];

(* Returns key->val pairs giving info about an intrade contract. *)
conInfo[ids_List] := Check[
  Import[itcUrl<>cat@@Riffle[ids,"&id="],"XML"][[2,3,All,2]], 
  prn["\n\n<h1>conInfo: no intrade data; check back in 10 minutes</h1>\n"];
  Exit[1]]

(* Turn the messy nested XML into a flat list of key->val pairs. *)
cleanPriceXML[ XMLElement["contractInfo", conInfo_, 
 {XMLElement["symbol", {}, {sym_}], 
  XMLElement["orderBook", {}, 
    {XMLElement["bids",{},bidList_], XMLElement["offers",{},askList_]}]}] ] :=
  Join[conInfo, {"bid"->If[bidList=={}, "0", bidList[[1,2,1,2]]], 
                 "ask"->If[askList=={}, "100", askList[[1,2,1,2]]]}]
cleanPriceXML[ XMLElement["contractInfo", conInfo_, {}] ] := conInfo

(* Returns key->val pairs giving latest price info for a contract. *)
priceInfo[ids_List] := Check[
  cleanPriceXML /@ Import[itpUrl<>cat@@Riffle[ids, "&id="], "XML"][[2,3]],
  prn["\n\n<h1>priceInfo: no intrade data; check back in 10 minutes</h1>\n"];
  Exit[1]]
  
(* Fetch some numeric value from the stuff returned by a call to intrade. *)
fetch[key_, info_List] := If[MemberQ[info[[All,1]], key],
  With[{s = key/.info}, If[SyntaxQ@s, eval@s, s]],
  Null]
     
(* Extract last trade price from price info, clipped to be betw bid & ask. *)
ltp[x_List] := With[{p = fetch["lstTrdPrc", x]},
  Which[
    p===Null, fetch["expiryPrice",x], 
    p==="-",  shoop[fetch["conID", x]], 
    True,     Clip[p, {fetch["bid",x], fetch["ask",x]}]]]
  

(*************************************************************************)
(*                Get Intrade Data and Run Simulations                   *)
(*************************************************************************)
  
(* Price history data (daily) for democrat/GOP presidential win. *)
datad = (phist /@ itd);
datar = (phist /@ itr);
  
(* Price history data for electoral vote count markets. *)
datadEC = (echist /@ DemECID);
datarEC = (echist /@ RepECID);

(* (no longer) data error: CA price on 9.11=3 instead of 91.2 *)
(* datad[[5,664,2;;-1]]=datad[[5,663,2;;-1]] *) 

(* data error from 2008: Sep 15 R price in MI is 98 *)
(*datar[[23, 657, 2;;-1]] = datar[[23, 656, 2;;-1]]; *)

(* List of daily closing prices. *)
closed = datad[[All,All,5]];
closer = datar[[All,All,5]];

returnsd = Table[Log[closed[[i,j-1]]/closed[[i,j]]], {i,51}, 
                 {j, -Min[WIND,Length[closed[1]]], -1}];
returnsr = Table[Log[closer[[i,j-1]]/closer[[i,j]]], {i,51}, 
                 {j, -Min[WIND,Length[closer[1]]], -1}];

(* days till the election *)
timeleft = Max[1,DateDifference[DateList[][[1;;3]], EDAY]];

(*
liquidityd = Table[Count[returnsd[[i, -WIND;;-1]], n_ /; n!=0], {i, 51}];
liquidityr = Table[Count[returnsr[[i, -WIND;;-1]], n_ /; n!=0], {i, 51}];
*)

(* pairwise state-by-state covariances, for democrat and republican prices *)
cvard = Table[Covariance[returnsd[[i]], returnsd[[j]]], {i,1,51}, {j,1,51}];
cvarr = Table[Covariance[returnsr[[i]], returnsr[[j]]], {i,1,51}, {j,1,51}];

closed = Take[#,-Min[WIND,Length[#]]]& /@ closed;
closer = Take[#,-Min[WIND,Length[#]]]& /@ closer;
  
nowd = ltp /@ priceInfo[itd];
nowr = ltp /@ priceInfo[itr];
  
(*  we're now getting the actual live prices instead of the last close 
nowd = closed[[All, -1]];
nowr = closer[[All, -1]];
*)

mleev = (Sign[nowd-nowr]/.{-1->0,0->1/2}).ev;
  
(* Map price (in [0,1]) to probability. *)
pp[p_] := Which[p<PMIN, 0, p>PMAX, 1, True, p]
  
(* Return a list of 51 D-win probabilities for a single simulation. *)
sample[] := Module[{ds, rs, endd, endr, chanced}, 
  (*
  ds = Exp[Transpose[RandomReal[mndist[Table[0,{51}], cvard], timeleft]]];
  rs = Exp[Transpose[RandomReal[mndist[Table[0,{51}], cvarr], timeleft]]];
  endd = nowd*Times @@@ ds;
  endr = nowr*Times @@@ rs;
  *)
  endd = nowd;
  endr = nowr;
  endd/(endd+endr)]

(* number of simulated outcomes for each poss number of electoral votes *)
counts = Table[0, {539}]; 
sbys = Table[0, {51}]; (* state-by-state: # of sims with D>.5 in the state *)
Do[s = sample[];
   i = 1+Sum[Boole[RandomReal[] < pp[s[[i]]]] * ev[[i]], {i, 51}];
   sbys += (Sign[s-1/2]/.{-1->0,0->1/2});
   counts[[i]]++, {NSIM}];
cdf = 1 - Accumulate@counts/Total@counts; (* 1 minus cdf of dem electrl votes *)
cdfr = 1 - Accumulate@Reverse@counts/Total@counts;
F[x_] := 1-cdf[[x+1]];  (* the real cdf -- eg, cdf[[1]] is Pr(D>=0) *)
Fr[x_] := 1-cdfr[[x+1]];
dwin = 1-F[269];  (* ie, Pr(D>=270) *)
rwin = 1-Fr[269];
  
(* Repeating the simulations but under pure independence... *)
countsi = Table[0, {539}];
Do[s = nowd/(nowd+nowr);
   i = 1+Sum[Boole[RandomReal[] < s[[i]]] * ev[[i]], {i, 51}];
   countsi[[i]]++, {NSIM}];
cdfi = 1-Accumulate@countsi/Total@countsi; (* 1 minus cdf of dem elctrl votes *)
cdfri = 1-Accumulate@Reverse@countsi/Total@countsi;
Fi[x_] := 1-cdfi[[x+1]];  (* the real cdf -- eg, cdf[[1]] is Pr(D>=0) *)
Fri[x_] := 1-cdfri[[x+1]];
dwini = 1-Fi[269];  (* ie, Pr(D>=270) *)
rwini = 1-Fri[269];

(* should keep things in terms of the more compact counts, not sims *)
sims = Flatten@MapIndexed[Table[#2[[1]]-1, {#1}]&, counts];
simdEC = Table[Count[sims, n_/;n>i]/Length[sims], {i,209,399,10}] * 100.
simrEC = Table[Count[538-sims, n_/;n>i]/Length[sims], {i,179,379,10}] * 100.

(* Show Number (get rid of the pesky trailing decimal point if any) *)
shn[x_] := StringReplace[ToString[If[IntegerQ@x,x,N@x]],re@"\\.$"->""]
  
indeplot[xmin_, xmax_] :=
  ListPlot[Select[Transpose[{Range[0,538], cdfi}], xmin<=#[[1]]-1<=xmax&], 
           Joined->False, PlotStyle->{LightGray}]
indeplotr[xmin_, xmax_] :=
  ListPlot[Select[Transpose[{Range[0,538], cdfri}], xmin<=#[[1]]-1<=xmax&], 
           Joined->False, PlotStyle->{LightGray}]
    

(* min and max for plot of cdf under independence of states *)
{ymini, ymaxi} = {Min[dwini,.01], Max[1-dwini,.99]};
xmini = With[{p=Position[1-cdfi, n_/;n>=ymini]}, 
             If[p=={}, 0, Min[269,p[[1,1]]-1]]];
xmaxi = With[{p=Position[1-cdfi, n_/;n<=ymaxi]},
             If[p=={}, 538, p[[-1,1]]-1]];
(* min and max for plot of cdf based on simulations *)             
{ymin, ymax} = {Min[dwin,.01], Max[1-dwin,.99]};
xmin = With[{p=Position[1-cdf, n_/;n>=ymin]}, If[p=={},0,Min[269,p[[1,1]]-1]]];
xmax = With[{p=Position[1-cdf, n_/;n<=ymax]}, If[p=={}, 538, p[[-1,1]]-1]];
{xmin, xmax} = {Min[xmin, xmini, 210], Max[xmax, xmaxi, 400]};
plotdata = Select[Transpose[{Range[0,538], cdf}], xmin<=#[[1]]-1<=xmax&];
leftdata = Select[plotdata, #[[1]]<269&];
tiedata = Select[plotdata, #[[1]]==269&];
rightdata = Select[plotdata, #[[1]]>=269&];
leftplot = If[leftdata=={}, Graphics[], 
  ListPlot[leftdata, Axes->False, PlotStyle->{Red, PointSize[ptsz]}]];
tieplot = If[tiedata=={}, Graphics[],
  ListPlot[tiedata, Axes->False, PlotStyle->{Purple, PointSize[ptsz]}]];
rightplot = If[rightdata=={}, Graphics[],
  ListPlot[rightdata, Axes->False, PlotStyle->{Blue, PointSize[ptsz]}]];
plot = Show[
  Graphics[{Green, Dashed, Line[{{269, 0}, {269, 1}}]}], 
  Plot[dwin, {x, xmin, xmax}, PlotStyle->{Green, Dashed}],
  indeplot[xmin, xmax], leftplot, rightplot, tieplot, 
  ListPlot[Transpose[{Range[210,400,10], datadEC[[All,-1,5]]/100}], 
           Filling->0],
  Axes->False, Frame->True, AspectRatio->1/GoldenRatio,
  (*PlotLabel->cat[
    "Probability of Obama (D) Getting At Least This Many Electoral Votes."],*)
  FrameLabel->{"Obama (D) Gets at Least This Many Electoral Votes", 
               "Probability"},
  BaseStyle->{FontFamily->"Geneva", FontSize->9}];
Export[path<>"cdf.png", plot, ImageSize->800];
system["scp ",path,"cdf.png "<>xpath];
  
  

(* repeat from mccain's perspective *)

(* min and max for plot of cdf under independence of states *)
{ymini, ymaxi} = {Min[rwini,.01], Max[1-rwini,.99]};
xmini = With[{p=Position[1-cdfri, n_/;n>=ymini]}, 
             If[p=={}, 0, Min[269,p[[1,1]]-1]]];
xmaxi = With[{p=Position[1-cdfri, n_/;n<=ymaxi]},
             If[p=={}, 538, p[[-1,1]]-1]];
(* min and max for plot of cdf based on simulations *)             
{ymin, ymax} = {Min[rwin,.01], Max[1-rwin,.99]};
xmin = With[{p=Position[1-cdfr, n_/;n>=ymin]}, 
             If[p=={},20,Min[269,p[[1,1]]-1]]];
xmax = With[{p=Position[1-cdfr, n_/;n<=ymax]}, 
             If[p=={}, 538, p[[-1,1]]-1]];
{xmin, xmax} = {Min[xmini, xmin, 180], Max[xmaxi, xmax, 380]};
{xmin, xmax} = {100, 380};
plotdata = Select[Transpose[{Range[0,538], cdfr}], xmin<=#[[1]]-1<=xmax&];
leftdata = Select[plotdata, #[[1]]<269&];
tiedata = Select[plotdata, #[[1]]==269&];
rightdata = Select[plotdata, #[[1]]>=269&];
leftplot = If[leftdata=={}, Graphics[], 
  ListPlot[leftdata, Axes->False, PlotStyle->{Blue, PointSize[ptsz]}]];
tieplot = If[tiedata=={}, Graphics[],
  ListPlot[tiedata, Axes->False, PlotStyle->{Purple, PointSize[ptsz]}]];
rightplot = If[rightdata=={}, Graphics[],
  ListPlot[rightdata, Axes->False, PlotStyle->{Red, PointSize[ptsz]}]];
plot = Show[
  Graphics[{Green, Dashed, Line[{{269, 0}, {269, 1}}]}], 
  Plot[rwin, {x, xmin, xmax}, PlotStyle->{Green, Dashed}],
  indeplotr[xmin, xmax], leftplot, rightplot, tieplot, 
  ListPlot[Transpose[{Range[180,380,10], datarEC[[All,-1,5]]/100}], 
           Filling->0],
  Axes->False, Frame->True, AspectRatio->1/GoldenRatio,
  FrameLabel->{"Romney (R) Gets at Least This Many Electoral Votes", 
               "Probability"},
  BaseStyle->{FontFamily->"Geneva", FontSize->9}];
Export[path<>"cdfr.png", plot, ImageSize->800];
system["scp ",path,"cdfr.png "<>xpath];
  



Export[path<>"pp.png", 
  Show[Plot[pp@x, {x,0,1}, PlotStyle->Green],
       Plot[x, {x,0,1}, PlotStyle->Dashed],
       pdata = Sort@Transpose[{nowd/(nowd+nowr), sbys/NSIM}];
       ListPlot[pdata, PlotStyle->PointSize@.015],
       ListPlot[pdata, Joined->True],
       Frame->True]];
system["scp ",path,"pp.png "<>xpath];

Export[path<>"bar1.png", 
  BarChart[Transpose@{datadEC[[All,-1,5]], simdEC}, 
           ChartLabels->{Table[i, {i,210,400,10}],None}],
  ImageSize->800];
system["scp ",path,"bar1.png "<>xpath];

Export[path<>"bar2.png", 
  BarChart[Transpose@{datarEC[[All,-1,5]], simrEC}, 
           ChartLabels->{Table[i, {i,180,380,10}],None}],
  ImageSize->800];
system["scp ",path,"bar2.png "<>xpath];
  
dneeds[_] = False;
rneeds[_] = False;
color[s_String] := Which[dneeds@s && rneeds@s, "purple", 
                         dneeds@s, "blue", rneeds@s, "red", True, "black"]
color[p_] := Which[p<40, "red", p>60, "blue", True, "purple"]

Module[{sum=0},
  Scan[(sum += #[[2]]; dneeds[#[[1]]] = True; If[sum>=269, Return[]])&,
       SortBy[Transpose[{states, ev, nowd/(nowd+nowr)}], -#[[3]]&]]; 
  sum=0;
  Scan[(sum += #[[2]]; rneeds[#[[1]]] = True; If[sum>=270, Return[]])&,
       SortBy[Transpose[{states, ev, nowr/(nowd+nowr)}], -#[[3]]&]]];

       
(* GENERATE HTML *)

str = OpenWrite[path<>"index.html"];
out = WriteString[str, ##, "\n"]&;

out["<html><head><title>Live Electoral Analysis</title></head><body>"];
out["<h1>Live Electoral Analysis (",
  DateString[{"Year","-","Month","-","Day"," ","Hour",":","Minute"," PDT"}],
  ") @T-",timeleft,"d</h1>"];
out["<h2>Cody Stumpo and Daniel Reeves</h2>"];
out["<img src=\"cdf.png\"><p>"];
out["The red/blue points give the probability distribution for number of ",
    "Obama electoral votes."];
out["The green dashed crosshairs give the probability of an Obama win."];
out["The bar chart gives the current probabilities ",
    "according to the electoral vote count markets on Intrade."];
out["The gray points give the probability distribution under the assumption ",
    "of perfect independence between states."];
out["<p><hr><p>"];
out["<img src=\"cdfr.png\"><p>"];
out["Same as above but from Romney's perspective."];
out["<p><hr><p>"];
out["<pre>"];
{mccain, obama} = (ltp /@ priceInfo[{743475, 743474}]);
out["Current Romney and Obama <a href=\"http://intrade.com\">intrade</a> prices: ", shn[mccain]," & ",shn[obama], 
    If[mccain+obama!=100, cat[" (renormalized: ", 
                             shn@Round[mccain/(mccain+obama)100,.1], " & ", 
                             shn@Round[obama/(mccain+obama)100,.1], ")"], ""]];
out[""];
out["Estimates of Obama's win probability:"];
out["  Based on simulations (and including the ",shn@N[F@269-F@268],
    " probability of a 269-269 tie): ", shn@N[1-F@268]];
out["  Assuming perfect correlation between states: ",
  Module[{sum=0},
    Scan[(sum += #[[1]]; If[sum>=269, Return[#[[2]]]])&, 
         SortBy[Transpose[{ev, nowd/(nowd+nowr)}], -#[[2]]&]]]];
out["  Assuming perfect independence between states: ", shn@N[1-Fi@268]];
out["Estimates of the number of Romney & Obama electoral votes:"];
out["  Expectation (mean of simulated distribution): ",
  shn@N[538 - counts.Range[0,538]/Total@counts], " & ",
  shn@N[counts.Range[0,538]/Total@counts]];
out["  Most frequently simulated outcome (mode of simulated distribution): ", 
  538 - Position[counts, Max@counts][[1,1]]-1, " & ",
  Position[counts, Max@counts][[1,1]]-1];
out["  Max likelihood estimate (\"leaning\" towards Romney & Obama): ", 
  shn@N[538 - mleev], " & ",
  shn@N@mleev];
out["  Expectation under the assumption of perfect independence betw states: ",
  shn@N[538 - countsi.Range[0,538]/Total@countsi], " & ",
  N[countsi.Range[0,538]/Total@countsi]];
out["</pre>"];
out["<hr>"];
out["<h3>Methodology</h3>"];
out["<ol>"];
out["<li>"];
out["Determine covariance matrix of 51 pairs of state electoral markets on"];
out["intrade, based on 90-day window of daily price movements.</li>"];
out["<li>Simulate Gaussian Copula of price movements per day, for as many"];
out["days as remain to the election.</li>"];
out["<li>Award each state's electoral votes according to ending price,"];
out["assuming limited fallibility of intrade prediction and independence of"];
out["prediction error:<br>price &lt; ",PMIN," : defeat,<br>",PMIN,
    " &lt; price &lt; ",PMAX," : victory assigned"];
out["with probability = price,<br>price &gt; ",PMAX," : victory.</li>"];
out["<li>Repeat simulation ",NSIM," times and gather statistics.</li>"];
out["</ol>"];
out["<hr>"];
out["The following table shows electoral votes by state, current intrade ",
    "prices, and percentage of simulations in which Obama won each state."];
out["States colored red (or purple) are the most likely set of states ",
    "sufficient for a Romney win."];
out["Similarly for states colored blue.  Probabilities/prices are colored ",
    "according to which candidate they favor."];
out["\n<p>"];
out["<table border=\"0\" cellspacing=\"0\" cellpadding=\"5\">"];
out["<tr><td>(Easiest way for Romney to win is to win these states...)</td>"];
out["<td>(Easiest way for Obama to win is to win these states...)</td></tr>"];
out["<tr valign=\"top\"><td width=\"50%\">"];
out["<table border=\"1\">"];   (* LEFT TABLE *)
out["<tr><th></th><th>EV</th><th>D-Price</th><th>R-Price</th>"];
out["<th>D+R</th><th>R/(D+R)</th><th>% R win</th></tr>"];
each[{s_, ev_, d_, r_, p_}, 
     SortBy[Select[
       Transpose[{states, ev, nowd, nowr, N@Round[sbys/NSIM*100,.1]}],
       rneeds[#[[1]]]&],
               #[[4]]/(#[[3]]+#[[4]])&],
  out["<tr><td><font color=",color[s],">",s,"</font></td><td>",ev,"</td>"];
  out["<td>",d,"</td><td>",r,"</td><td>",d+r,"</td><td><font color=",
      color[d/(d+r)*100],">",shn@N@Round[r/(d+r)*100,.1],
      "</font></td><td><font color=",color[p],">",100-p,"</font></td></tr>"]];
out["</table>"];  (* END LEFT TABLE *)
out["</td><td width=\"50%\">"];
out["<table border=\"1\">"];   (* RIGHT TABLE *)
out["<tr><th></th><th>EV</th><th>D-Price</th><th>R-Price</th>"];
out["<th>D+R</th><th>D/(D+R)</th><th>% D win</th></tr>"];
each[{s_, ev_, d_, r_, p_}, 
     SortBy[Select[
       Transpose[{states, ev, nowd, nowr, N@Round[sbys/NSIM*100,.1]}],
       dneeds[#[[1]]]&],
               #[[3]]/(#[[3]]+#[[4]])&],
  out["<tr><td><font color=",color[s],">",s,"</font></td><td>",ev,"</td>"];
  out["<td>",d,"</td><td>",r,"</td><td>",d+r,"</td><td><font color=",
      color[d/(d+r)*100],">",shn@N@Round[d/(d+r)*100,.1],
      "</font></td><td><font color=",color[p],">",p,"</font></td></tr>"]];
out["</table>"];  (* END RIGHT TABLE *)
out["</td></tr></table>"];
out["<p><hr><p>"];
out["<img src=\"pp.png\"><p>"];
out["Given current intrade price (x axis) and the price history (across 51 ",
    "electoral markets), what is the chance (y axis) that the price on ",
    "election day will be above 50?"];
out["<p><hr><p>"];
(*
out["<img src=\"indep.png\">"];
out["<p><hr><p>"];
*)
out["<img src=\"bar1.png\"><br>"];
out["What is the simulated chance Obama will win x or more electoral
votes (green) and what is the Intrade contract price on the same
question (blue)"];
out["<p><hr><p>"];
out["<img src=\"bar2.png\"><br>"];
out["What is the simulated chance Romney will win x or more electoral ",
    "votes (green) and what is the Intrade contract price on the same ",
    "question (blue)"];
out["<p><hr>This page was generated in Mathematica at ",
  DateString[{"Year","-","Month","-","Day"," ","Hour",":","Minute"," PDT"}]];
  
Close@str;
system["scp ",path,"index.html "<>xpath];


(***************************** SCRATCH **************************************)
  
(*
prn["<h3>Scratch Area -- Ignore Everything Below</h3>"];
prn[N@SortBy[Transpose[{states,ev,nowd/(nowd+nowr)}],Abs[#[[3]]-1/2]&]];
prn["</body></html>"]; 
*)
  
(* translate back to a list of simulated outcomes *)
(* sims = Flatten@MapIndexed[Table[#2[[1]] - 1, {#1}] &, counts]; *)
 
(*
Needs["Histograms`"];
hr = Histogram[Select[sims, # < 269 &], 
       Ticks -> {Table[10 x, {x, 21, 41}], None}, 
       HistogramCategories -> Table[x, {x, 538}], BarStyle->Red];
hd = Histogram[Select[sims, #>268&], 
       Ticks -> {Table[10 x, {x, 21, 41}], None}, 
       HistogramCategories -> Table[x, {x, 538}], BarStyle->Blue];
hist = Show[hr, hd, PlotRange -> {{220, 400}, {0, Max[counts]}}];

chanced = Table[closed[[i, -WIND;;-1]]/(closed[[i, -WIND;;-1]] + 
                                          closer[[i, -WIND;;-1]]), {i, 51}];

swingevD = 
  Total[Extract[ev, Position[chanced[[All,-1]], n_ /; .5<n<.7]]];

ohpa = Table[Switch[Total[simsdet[[i, {36, 39}]]], 
               0, "Neither", 20, "OH", 21, "PA", 41, "Both"], {i, 5000}];
*)

(* Very little chance of winning Ohio but not Pennsylvania. 
   Probabilities very different from if OH & PA were independent. *)

(*
chanced = 
  Table[closed[[i, -WIND;;-1]]/(closed[[i, -WIND;;-1]] + 
                                  closer[[i, -WIND;;-1]]), {i, 51}];

sbys = {Table[Count[simsdet[[All,i]], n_ /; n>0]/Length[simsdet], {i, 51}]//N, 
        detailIDdr[[1]], chanced[[All,-1]]//N}; 
*)

  
(*
returnsc = 
  Table[Log[chanced[[i, j-1]]/chanced[[i, j]]], {i, 51}, {j, -WIND, -1}];
winning = Table[Total[Table[If[chanced[[i,-j]]>.5, ev[[i]], 0], {i, 51}]], 
                {j, WIND}];
DateListPlot[Reverse[winning], {2008,7,4}, Joined->True, 
                                           PlotRange->{Automatic, {250, 350}}]
*)


(*
chanced = Table[closed[[i, -WIND;;-1]]/(closed[[i, -WIND;;-1]] + 
                                          closer[[i, -WIND;;-1]]), {i, 51}]; 
     
today = Total[Table[If[chanced[[i,-1]]>.5, ev[[i]], 0], {i, 51}]];
*)

(* Everything breaking exactly like today is the mode of the distribution *)
(* counts[[today]] == mode; *)
