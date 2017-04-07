{
//*-*-*-*-*-*-*Set visibility for this node and its sons*-*-*-*-*--*-*-*-*-*-*
//*-*          =========================================
//*-*  vis = 3  node is drawn and its sons are drawn
//*-*  vis = 2  node is not drawn but its sons are drawn
//*-*  vis = 1  (default) node is drawn
//*-*  vis = 0  node is not drawn
//*-*  vis = -1 node is not drawn. Its sons are not drawn
//*-*  vis = -2 node is drawn. Its sons are not drawn
//*-*  vis = -3 Only node leaves are drawn
//*-*  vis = -4 Node is not drawn. Its immediate sons are drawn
//*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

HALL1->ImportShapeAttributes();
HALL1->SetLineColor(0);
//  Set Node attributes
HALL1->SetVisibility(2);   //node is not drawn but its sons are drawn
ACR11->SetVisibility(2);
ACS21->SetVisibility(2); 
ACS31->SetVisibility(2); 
AHM11->SetVisibility(2);
ASP21->SetVisibility(2);
ASPM1->SetVisibility(2);
//MAG21->SetVisibility(2);
MTAR1->SetVisibility(2);
SDC_1->SetVisibility(2);
SDC_2->SetVisibility(2);
SM1M1->SetVisibility(2);
SMM_1->SetVisibility(2);
SMM_2->SetVisibility(2);
SMM_3->SetVisibility(2);
VSB_1->SetVisibility(2);
VSB_2->SetVisibility(2);


MAG21->SetVisibility(1);//SM2
MAG21->SetLineColor(0);
M2GP1->SetLineColor(0);
M2D11->SetLineColor(0);
M2DH1->SetLineColor(0);

//-------PT---------------------------------
PTSU1->SetVisibility(2); //Polarized Target

TSO11->SetVisibility(2);
TSO21->SetVisibility(2);
TSOL1->SetVisibility(2);

PTSX1->SetVisibility(2);

CCON1->SetVisibility(2);
CTUB1->SetVisibility(2);
FIBG1->SetVisibility(2);
FLAA1->SetVisibility(2);
FLAF1->SetVisibility(2);
FLAS1->SetVisibility(2);
FOCU1->SetVisibility(2);
FOI11->SetVisibility(2);
FOI21->SetVisibility(2);
FOI31->SetVisibility(2);
FOI41->SetVisibility(2);
FOI51->SetVisibility(2);
FOI61->SetVisibility(2);
OBL11->SetVisibility(2);
OBL21->SetVisibility(2);
HELI1->SetVisibility(2);
HNKM1->SetVisibility(2);

UPST1->SetVisibility(-2);//Target Material
DWST1->SetVisibility(-2);//Target Material
UPST1->SetLineColor(2);
DWST1->SetLineColor(2);
//----------------------------------

//-----------SM1--------------------
SM1M1->SetLineColor(0);
FLDM1->SetLineColor(0);
POLT1->SetLineColor(0);
POLB1->SetLineColor(0);
//----------------------------------

//-----------RICH-------------------
RCH11->SetVisibility(-2); //RICH1
RCH11->SetLineColor(4);

//-----------MUON FILTER------------
MU1F1->SetVisibility(-2);
MU2F1->SetVisibility(-2);
MU1F1->SetLineColor(7);
MU2F1->SetLineColor(7);
//----------------------------------

//----------------------------------
ACS21->SetVisibility(-1); //Not Draw
ACS31->SetVisibility(-1); 
SMM_1->SetVisibility(-1);
SMM_2->SetVisibility(-1);
SMM_3->SetVisibility(-1);
SDC_1->SetVisibility(-1);
SDC_2->SetVisibility(-1);
VSB_1->SetVisibility(-1);
VSB_2->SetVisibility(-1);

MG1_5->SetVisibility(-1);
MG1_6->SetVisibility(-1);
MG1_7->SetVisibility(-1);
MG1_8->SetVisibility(-1);
MG1_9->SetVisibility(-1);

//------- OMEGA A-Chamber ---------- 
AAA_1->SetVisibility(-2);
AAA_2->SetVisibility(-2);
AAA_3->SetVisibility(-2);
AAA_4->SetVisibility(-2);
AAA_5->SetVisibility(-2);
AAA_6->SetVisibility(-2);
AAA_1->SetLineColor(0);
AAA_2->SetLineColor(0);
AAA_3->SetLineColor(0);
AAA_4->SetLineColor(0);
AAA_5->SetLineColor(0);
AAA_6->SetLineColor(0);
//----------------------------------

//------- OMEGA B-Chamber ---------- 
BB1_1->SetVisibility(-2);
BB1_2->SetVisibility(-2);
BB2_1->SetVisibility(-2);
BB2_2->SetVisibility(-2);
BB1_1->SetLineColor(0);
BB1_2->SetLineColor(0);
BB2_1->SetLineColor(0);
BB2_2->SetLineColor(0);
//-----------------------------------

//------- Hodoscope -----------------
HOD11->SetVisibility(-1);
HOD21->SetVisibility(-1);
HOD31->SetVisibility(-2);
HOD41->SetVisibility(-2);
HOD51->SetVisibility(-2);
HOD31->SetLineColor(5);
HOD41->SetLineColor(5);
HOD51->SetLineColor(5);

H4VS1->SetVisibility(-2);
H4VV1->SetVisibility(-2);
H5VS1->SetVisibility(-2);
H4VS1->SetLineColor(5);
H4VV1->SetLineColor(5);
H5VS1->SetLineColor(5);

//HO1H1->SetVisibility(0);
//HO2H1->SetVisibility(0);
//-----------------------------------

VSD_1->SetVisibility(-1);
VSD_2->SetVisibility(-1);
VSD_3->SetVisibility(-1);
VSD_4->SetVisibility(-1);
//VSB_2->SetVisibility(0);



//-------Muon Wall ------------------
//MU1W1->SetVisibility(-2);
MU2W1->SetVisibility(-2);
MU3W1->SetVisibility(-2);
//MU1W1->SetLineColor(0);
MU2W1->SetLineColor(11);
MU3W1->SetLineColor(11);
//-----------------------------------



//-------Hadron Calorimeter-------------------
HC1_1->SetVisibility(-2);
HC2_1->SetVisibility(-2);
HC1_1->SetLineColor(6);
HC2_1->SetLineColor(6);
//--------------------------------------------

//ASP21->SetVisibility(-1);

}

// void SetLineColor(TNode *node,Int_t color){
//   node->SetLineColor(color);
//   TList *nodeList = node->GetListOfNodes();
//   TList *subNodeList;
//   TNode *subNode1,*subNode2;
//   TIter nNode(nodeList);
//   while((subNode1 = (TNode*)nNode())) {
//     subNode1->SetLineColor(0);
//     if((subNodeList = subNode1->GetListOfNodes())){
//       TIter nSubNode(subNodeList);
//       while((subNode2 = (TNode*)nSubNode())) {
// 	SetLineColor(subNode2,color);
//       }          
//     }
//   }
// }



