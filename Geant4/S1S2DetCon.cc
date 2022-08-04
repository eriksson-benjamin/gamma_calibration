// Anders Hjalmarsson
//
// Geometry for the S1S2 detector 1999/12/02 09.15
// GEANT4 program
//
// $Id: S1S2DetCon.cc

#include "S1S2DetCon.hh"

#include "S1S2DetMessenger.hh"

#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "globals.hh"

#include "G4MaterialTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpBoundaryProcess.hh"

//#include "S1S2SD.hh"

#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

S1S2DetCon::S1S2DetCon():
solidWorld(NULL), logicalWorld(NULL), physiWorld(NULL), solidS1(NULL),
logicalS1(NULL), physiS1(NULL), solidS2(NULL), logicalS2(NULL), physiS2(NULL),
BC404(NULL),Air(NULL),boolLightGuide("on")
{
  x1Trd = 2.5*cm;
  x2Trd = 2.5*cm;
  zTrd = 35*cm;

  tiltangle = 0.0*deg;
  
  meanSangle = 27.5*deg;
  rangle = (180*deg-2*meanSangle);
  meanFlightpath = sqrt(2*pow(704.6,2)*(1-cos(rangle)))*mm;
  
  xposS2Abs = meanFlightpath*cos(meanSangle);
  yposS2Abs = meanFlightpath*sin(meanSangle);
  zposS2Abs = meanFlightpath*sin(meanSangle);
  
  hightS1 = 1.25*cm;

  outRadS1 = 1.8*cm;
  
  NumberOfS1 = 5;
  NumberOfS2 = 40;
  
  //y1Trd = 2*(yposS2Abs-zTrd/2*cos(angle1))*3.1416/NumberOfS2-0.5*cm;
  //y2Trd = 2*(yposS2Abs+zTrd/2*cos(angle1))*3.1416/NumberOfS2-0.5*cm;

  y1Trd = 95.0*mm;
  y2Trd = 134.7*mm;
  detMessenger = new S1S2DetMessenger(this);
}

S1S2DetCon::~S1S2DetCon()
{
  delete detMessenger;
}

G4VPhysicalVolume* S1S2DetCon::Construct()
{
  DefineMaterials();
  return ConstructAll();
}

G4VPhysicalVolume* S1S2DetCon::ConstructAll()
{
  //---- World -----------------------------------------
  xWorld = 400*cm;
  G4double yWorld = 400*cm;
  G4double zWorld = yWorld;
  
  solidWorld = new G4Box("World", 0.5*xWorld, 0.5*yWorld, 0.5*zWorld);
  logicalWorld = new G4LogicalVolume(solidWorld, Air, "WorldLV", 0, 0, 0);
  physiWorld = new G4PVPlacement(0, G4ThreeVector(), "WorldPV",
				 logicalWorld, NULL, false, 0);
  
  //---- S1 detector consists of 5 scintillators
  solidS1 = NULL, logicalS1 = NULL, physiS1 = NULL;
  G4double inRadS1 = 0.*cm;
  G4double startAngS1 = 0.*deg;
  G4double spanAngS1 = 360.*deg;
  
  solidS1 = new G4Tubs("S1", inRadS1, outRadS1, 0.5*hightS1, startAngS1, spanAngS1);
  logicalS1 = new G4LogicalVolume(solidS1, BC404, "S1LV", 0, 0, 0);

 
  G4RotationMatrix RotTub;
  RotTub.rotateY(90.*deg);
  
  
  RotTub.rotateY(90.*deg);
  
  for (G4int i = 0; i < NumberOfS1; i++)
    //for (G4int i = 0; i < 1; i++)
    {
      G4double xposS1;
      G4double yposS1 = 0;
      G4double zposS1 = 0;
      G4double angle = 90.*deg;
      xposS1 = (hightS1/cm+0.1)*(i-2)*cm;
      //xposS1 = hightS1*(i-0.5);
      //xposS1 = 0;
      G4RotationMatrix RotMat;
      RotMat.rotateY(angle);
      
      physiS1 = new G4PVPlacement(G4Transform3D(RotMat, 
						G4ThreeVector(xposS1, yposS1, zposS1)),
				  "S1PV", logicalS1, physiWorld, false, i);
    }
  
  //---- S2 detector consists of NumberOfS2 scintillators ----------------
  //---- Size 35 * 10.0/6.5 * 1.25 cm and area 12.5 cm^2 to the PM-tube -- 
  solidS2 = NULL, logicalS2 = NULL, physiS2 = NULL;
  solidS2 =new G4Trd("S2", 0.5*x1Trd, 0.5*x2Trd, 0.5*y1Trd, 0.5*y2Trd, 0.5*zTrd);
  logicalS2 = new G4LogicalVolume(solidS2, BC404, "S2LV", 0, 0, 0);
  
  
  solidLGBox=NULL, solidLGTubs=NULL, solidLGS1=NULL, logicalLGS1=NULL, physiLGS1=NULL;
  G4RotationMatrix RotLG;
  RotLG.rotateY(90.*deg);
  G4RotationMatrix RotLG_b;
  RotLG_b.rotateX(180.0*deg);
  G4RotationMatrix RotLG_c;
  RotLG_c.rotateZ(90.*deg);

  G4double LGwidth = outRadS1*sin(120*deg)/sin(30*deg);
  solidLGBox = new G4Box("LGBox", 0.5*5.0*mm, 0.5*LGwidth*mm, 0.5*84.0*mm);
  solidLGTubs = new G4Tubs("LGTubs", inRadS1, outRadS1, 0.51*hightS1, startAngS1, spanAngS1);
  
  solidLGTrd = new G4Trd("LGTrd", 0.5*5.0*mm, 0.5*15*mm, 0.5*LGwidth*mm, 0.5*15*mm, 0.5*25*mm);
  
  solidLGCons = new G4Cons("LGCons", 0.5*(LGwidth+0.1)*mm, 0.5*(LGwidth+10.0)*mm, 0.5*15*mm, 0.5*25*mm,0.5*25.1*mm,0*deg,360*deg);
  
  solidLGS1 = new G4SubtractionSolid("solidLGS1",solidLGBox,solidLGTubs,G4Transform3D(RotLG,G4ThreeVector(0, 0, 42*mm)));
  solidLGS1_b = new G4SubtractionSolid("solidLGS1_b",solidLGTrd,solidLGCons,0,G4ThreeVector(0,0,0));
  logicalLGS1 = new G4LogicalVolume(solidLGS1, LG, "LGS1LV", 0, 0, 0);
  logicalLGS1_b = new G4LogicalVolume(solidLGS1_b, LG, "LGS1LV_b", 0, 0, 0);

  solidLGTubs_b = new G4Tubs("LGTubs_b", 0, 0.5*15.0*mm, 0.5*50*mm, startAngS1, spanAngS1);
  logicalLGS1_c = new G4LogicalVolume(solidLGTubs_b, LG, "LGS1LV_c", 0, 0, 0);
  
  G4RotationMatrix RotSolidLG;
  G4double angleLG=0;

  G4int n=0;

  if(boolLightGuide == "on")
    {
      for(G4int j=0;j<NumberOfS1;j++)
	{
	  G4double posLGx=(hightS1/cm+0.1)*(j-2)*cm;
	  for(G4int i=0;i<3;i++)
	    {
	      G4double posLGy=42.0*sin(angleLG)*mm;
	      G4double posLGz=-42.0*cos(angleLG)*mm;
	      G4double posLGy_b=(84.0+12.5)*sin(angleLG)*mm;
	      G4double posLGz_b=(-84.0-12.5)*cos(angleLG)*mm;
	      G4double posLGy_c=(84.0+25.0+25.0)*sin(angleLG)*mm;
	      G4double posLGz_c=(-84.0-25.0-25.0)*cos(angleLG)*mm;
	      physiLGS1_b = new G4PVPlacement(G4Transform3D(RotLG_b,G4ThreeVector(posLGx,posLGy_b,posLGz_b)),
					      "LGS1PV_b",logicalLGS1_b,physiWorld,false,n);
	      physiLGS1 = new G4PVPlacement(G4Transform3D(RotSolidLG,G4ThreeVector(posLGx,posLGy,posLGz)),
					    "LGS1PV",logicalLGS1,physiWorld,false,n);
	      physiLGS1_c = new G4PVPlacement(G4Transform3D(RotLG_c,G4ThreeVector(posLGx,posLGy_c,posLGz_c)),
					      "LGS1PV_c",logicalLGS1_c,physiWorld,false,n);
	      angleLG += 120*deg;
	      RotSolidLG.rotateX(120*deg);
	      RotLG_b.rotateX(120*deg);
	      RotLG_c.rotateX(120*deg);
	      n++;
	    }
	  G4double angleLG2 = 360*deg/(NumberOfS1*3);
	  angleLG += angleLG2;
	  RotSolidLG.rotateX(angleLG2);
	  RotLG_b.rotateX(angleLG2);
	  RotLG_c.rotateX(angleLG2);
	}
    }
  
  solidPMS1 = new G4Tubs("PMS1Tubs",0.5*32.72*mm,0.5*34.0*mm,0.5*151*mm,0*deg,360*deg);
  logicalPMS1 = new G4LogicalVolume(solidPMS1, PM,"PMS1LV",0,0,0);
  //physiPMS1 = new G4PVPlacement(0,G4ThreeVector(40*mm,0,0),"PMS1PV",logicalPMS1,physiWorld,false,0);
  
  G4double anglePM=0;
  G4RotationMatrix RotPM;
  RotPM.rotateX(180.0*deg);
  
  n=0;
  for(G4int j=0;j<NumberOfS1;j++)
    {
      G4double posPMx=(hightS1/cm+0.1)*(j-2)*cm;
      for(G4int i=0;i<3;i++)
	{
	  G4double posPMy=(84.0+25.0+50.0+0.5*137)*sin(anglePM)*mm;
	  G4double posPMz=(-84.0-25-50.0-0.5*137)*cos(anglePM)*mm;
	  physiPMS1 = new G4PVPlacement(G4Transform3D(RotPM,G4ThreeVector(posPMx,posPMy,posPMz)),
					"PMS1PV",logicalPMS1,physiWorld,false,n);
	  anglePM += 120*deg;
	  RotPM.rotateX(120*deg);
	}
      G4double anglePM2 = 360*deg/(NumberOfS1*3);
      anglePM += anglePM2;
      RotPM.rotateX(anglePM2);
    }
	  
  for (G4int i = 0; i < NumberOfS2; i++)
    {
      G4double angle2;
      angle2 = 360*deg*i/NumberOfS2;
      if(NumberOfS2 == 2)
	{
	  angle2 = 180*deg*i;
	}
      
      G4RotationMatrix RotMat;
      RotMat.rotateX(-90*deg);
      RotMat.rotateZ(angle1);
      RotMat.rotateX(angle2);
      G4double xposS2 = xposS2Abs;
      G4double yposS2 = yposS2Abs*cos(angle2);
      G4double zposS2 = zposS2Abs*sin(angle2);
      physiS2 = new G4PVPlacement(G4Transform3D(RotMat, G4ThreeVector(xposS2, yposS2, zposS2)),
				  "S2PV", logicalS2, physiWorld, false, i);

      
    }
  
  solidLGS2Trd = new G4Trd("LGS2Trd", 0.5*x1Trd, 0.5*50*mm,0.5*y2Trd , 0.5*50*mm, 0.5*100*mm);
  solidLGS2Cons = new G4Cons("LGS2Cons", 0.5*(y2Trd+5.0)*mm, 0.5*(y2Trd+25.0)*mm, 0.5*50*mm, 0.5*75*mm,0.5*101*mm,0*deg,360*deg);
  solidLGS2 = new G4SubtractionSolid("solidLGS2",solidLGS2Trd,solidLGS2Cons,0,G4ThreeVector(0,0,0));
  solidLGS2Tubs = new G4Tubs("LGS2Tubs",0,0.5*50.8*mm,0.5*50*mm,0*deg,360*deg);

  logicalLGS2 = new G4LogicalVolume(solidLGS2, LG, "LGS2LV", 0, 0, 0);
  logicalLGS2_b = new G4LogicalVolume(solidLGS2Tubs,LG,"LGS2LV_b",0,0,0);

  if(boolLightGuide == "on")
    {
      for(G4int i = 0; i < NumberOfS2; i++)
	{
	  G4double angle2;
	  
	  angle2 = 360*deg*i/NumberOfS2;
	  G4RotationMatrix RotLGS2;
	  RotLGS2.rotateX(-90*deg);
	  RotLGS2.rotateZ(angle1);
	  RotLGS2.rotateX(angle2);
	  G4double xposLG = xposS2Abs-((zTrd/2+50)*sin(angle1))*mm;
	  G4double yposLG = (yposS2Abs+((zTrd/2+50)*cos(angle1))*mm)*cos(angle2);
	  G4double zposLG = (zposS2Abs+((zTrd/2+50)*cos(angle1))*mm)*sin(angle2);
	  physiLGS2 = new G4PVPlacement(G4Transform3D(RotLGS2,G4ThreeVector(xposLG, yposLG, zposLG)),
					"LGS2PV",logicalLGS2,physiWorld, false, i);
	  G4RotationMatrix RotLGS2Tub;
	  RotLGS2Tub.rotateX(-90*deg);
	  RotLGS2Tub.rotateZ(angle1);
	  RotLGS2Tub.rotateX(angle2);
	  G4double xposLGTub = xposS2Abs-((zTrd/2+100+25)*sin(angle1))*mm;
	  G4double yposLGTub = (yposS2Abs+((zTrd/2+100+25)*cos(angle1))*mm)*cos(angle2);
	  G4double zposLGTub = (zposS2Abs+((zTrd/2+100+25)*cos(angle1))*mm)*sin(angle2);
	  physiLGS2_b = new G4PVPlacement(G4Transform3D(RotLGS2Tub,G4ThreeVector(xposLGTub, yposLGTub, zposLGTub)),
					  "LGS2PV_b",logicalLGS2_b,physiWorld, false, i);
	  
	}			   
    }
  
  solidPMS2 = new G4Tubs("PMS2Tubs",0.5*67.0*mm,0.5*70.0*mm,0.5*222.5*mm,0*deg,360*deg);
  logicalPMS2 = new G4LogicalVolume(solidPMS2, PM,"PMS2LV",0,0,0);
  
  solidPMS2SH = new G4Tubs("solidPMS2SH",0.5*82.6*mm,0.5*88.9*mm,0.5*224*mm,0*deg,360*deg);
  solidPMS2SHS = new G4Tubs("solidPMS2SHS",0.5*86.0*mm,0.5*90.0*mm,0.5*57.0*mm,0*deg,360*deg);
  logicalPMS2SH = new G4LogicalVolume(solidPMS2SH,PM,"solidPMS2SHLV",0,0,0);
  logicalPMS2SHS = new G4LogicalVolume(solidPMS2SHS,PM,"solidPMS2SHSLV",0,0,0);
  
  
  
  for(G4int i = 0; i < NumberOfS2; i++)
    {
      G4double angle2;
      angle2 = 360*deg*i/NumberOfS2;
      G4RotationMatrix RotPMS2;
      RotPMS2.rotateX(-90*deg);
      RotPMS2.rotateZ(angle1);
      RotPMS2.rotateX(angle2);
      G4double xposPM=xposS2Abs-((zTrd/2+100+50+192.5*0.5)*sin(angle1))*mm;
      G4double yposPM=(yposS2Abs+((zTrd/2+100+50+192.5*0.5)*cos(angle1))*mm)*cos(angle2);
      G4double zposPM= (zposS2Abs+((zTrd/2+100+50+192.5*0.5)*cos(angle1))*mm)*sin(angle2);
      physiPMS2 = new G4PVPlacement(G4Transform3D(RotPMS2,G4ThreeVector(xposPM, yposPM, zposPM)),
				    "PMS2PV",logicalPMS2,physiWorld, false, i);
      G4double xposSH = xposS2Abs-((zTrd/2+100+50+194*0.5)*sin(angle1))*mm;
      G4double yposSH = (yposS2Abs+((zTrd/2+100+50+194*0.5)*cos(angle1))*mm)*cos(angle2);
      G4double zposSH = (zposS2Abs+((zTrd/2+100+50+194*0.5)*cos(angle1))*mm)*sin(angle2);
      physiPMS2SH = new G4PVPlacement(G4Transform3D(RotPMS2,G4ThreeVector(xposSH,yposSH,zposSH)),
					    "solidPMS2SHPV",logicalPMS2SH,physiWorld, false, i);
      G4double xposSHS = xposS2Abs-((zTrd/2+100+50-30-8.5)*sin(angle1))*mm;
      G4double yposSHS = (yposS2Abs+((zTrd/2+100+50-30-8.5)*cos(angle1))*mm)*cos(angle2);
      G4double zposSHS = (zposS2Abs+((zTrd/2+100+50-30-8.5)*cos(angle1))*mm)*sin(angle2);
      physiPMS2SH = new G4PVPlacement(G4Transform3D(RotPMS2,G4ThreeVector(xposSHS,yposSHS,zposSHS)),
				      "solidPMS2SHSPV",logicalPMS2SHS,physiWorld, false, i);
    }

  G4RotationMatrix RotLGS1Ring;
  RotLGS1Ring.rotateY(90*deg);
  G4RotationMatrix RotPMS1Ring;
  RotPMS1Ring.rotateY(90*deg);
  G4RotationMatrix RotLGS1RingSub;
  RotLGS1RingSub.rotateY(90*deg);
  G4RotationMatrix RotPMS1RingSub;
  RotPMS1RingSub.rotateY(90*deg);
  G4double angleLGS1Ring=0*deg;
  

  solidLGS1Ring = new G4Tubs("LGS1Ring",120.0*mm, 130.0*mm,0.5*100.0*mm,0*deg,360*deg);
  solidLGS1Ring_b = new G4Tubs("LGS1Ring_b",0.0*mm,0.5*15.2*mm,0.5*15*mm,0.0*deg,360.0*deg);
  
  solidPMS1Ring_a = new G4Tubs("PMS1Ring_a",157.01*mm, 177.01*mm,0.5*50.0*mm,0*deg,360*deg);
  solidPMS1Ring_b = new G4Tubs("PMS1Ring_b",248.01*mm,268.01*mm,0.5*50.0*mm,0*deg,360*deg);
  solidPMS1Ring_c = new G4Tubs("PMS1Ring_c",0.0*mm, 0.5*24.0*mm,0.5*40.0*mm,0*deg,360*deg);

  n=0;
  for(G4int j=0;j<NumberOfS1;j++)
    {
      G4double posLGS1Ringz=(hightS1/cm+0.1)*(j-2)*cm;
      for(G4int i=0;i<3;i++)
	{
	  G4double posLGS1Ringy=125.0*sin(angleLGS1Ring)*mm;
	  G4double posLGS1Ringx=125.0*cos(angleLGS1Ring)*mm;
       
	  G4double posPMS1Ringy=167.01*sin(angleLGS1Ring)*mm;
	  G4double posPMS1Ringx=167.01*cos(angleLGS1Ring)*mm;
	  
	  G4double posPMS1Ring2y=258.01*sin(angleLGS1Ring)*mm;
	  G4double posPMS1Ring2x=258.01*cos(angleLGS1Ring)*mm;
	  
	  if(i==0 && j==0)
	    {
	      solidLGS1Ring_c[n] = new G4SubtractionSolid("solidLGS1Ring_c",solidLGS1Ring,solidLGS1Ring_b,
							  G4Transform3D(RotLGS1RingSub,G4ThreeVector(posLGS1Ringx,posLGS1Ringy,posLGS1Ringz)));
	      solidPMS1Ring_d[n] = new G4SubtractionSolid("solidPMS1Ring_d",solidPMS1Ring_a,solidPMS1Ring_c,
							   G4Transform3D(RotPMS1RingSub,G4ThreeVector(posPMS1Ringx,posPMS1Ringy,posLGS1Ringz)));
	      solidPMS1Ring_e[n] = new G4SubtractionSolid("solidPMS1Ring_e",solidPMS1Ring_b,solidPMS1Ring_c,
							  G4Transform3D(RotPMS1RingSub,G4ThreeVector(posPMS1Ring2x,posPMS1Ring2y,posLGS1Ringz)));
	    }
	  else
	    {
	      solidLGS1Ring_c[n]=new G4SubtractionSolid("solidLGS1Ring_c",solidLGS1Ring_c[n-1],solidLGS1Ring_b,
							  G4Transform3D(RotLGS1RingSub,G4ThreeVector(posLGS1Ringx,posLGS1Ringy,posLGS1Ringz)));
	      solidPMS1Ring_d[n] = new G4SubtractionSolid("solidPMS1Ring_d",solidPMS1Ring_d[n-1],solidPMS1Ring_c,
							   G4Transform3D(RotPMS1RingSub,G4ThreeVector(posPMS1Ringx,posPMS1Ringy,posLGS1Ringz)));
	      solidPMS1Ring_e[n] = new G4SubtractionSolid("solidPMS1Ring_e",solidPMS1Ring_e[n-1],solidPMS1Ring_c,
							   G4Transform3D(RotPMS1RingSub,G4ThreeVector(posPMS1Ring2x,posPMS1Ring2y,posLGS1Ringz)));
	    }
	  angleLGS1Ring += 120*deg;
	  RotLGS1RingSub.rotateZ(120*deg);
	  RotPMS1RingSub.rotateZ(120*deg);
	  n++;
	}
      G4double angleLGS1Ring2 = 360*deg/(NumberOfS1*3);
      angleLGS1Ring += angleLGS1Ring2;
      RotLGS1RingSub.rotateZ(angleLGS1Ring2);
      RotPMS1RingSub.rotateZ(angleLGS1Ring2);
    }
  
  logicalLGS1Ring = new G4LogicalVolume(solidLGS1Ring_c[n-1],Al,"LGS1RingLV",0,0,0);

  logicalPMS1Ring_a= new G4LogicalVolume(solidPMS1Ring_d[n-1],PM,"LGS1Ring_aLV",0,0,0);
  
  logicalPMS1Ring_b= new G4LogicalVolume(solidPMS1Ring_e[n-1],PM,"LGS1Ring_bLV",0,0,0);

  physiLGS1Ring = new G4PVPlacement(G4Transform3D(RotLGS1Ring,G4ThreeVector(0.0*mm,0.0*mm,0.0*mm)),
  				    "LGS1RingPV",logicalLGS1Ring,physiWorld, false, 0);
  
  //physiPMS1Ring_a = new G4PVPlacement(G4Transform3D(RotPMS1Ring,G4ThreeVector(0.0*mm,0.0*mm,0.0*mm)),
  //			      "PMS1Ring_aPV",logicalPMS1Ring_a,physiWorld, false, 0);
  
  //physiPMS1Ring_b = new G4PVPlacement(G4Transform3D(RotPMS1Ring,G4ThreeVector(0.0*mm,0.0*mm,0.0*mm)),
  //			      "PMS1Ring_bPV",logicalPMS1Ring_b,physiWorld, false, 0);
  
  
  
  G4RotationMatrix RotLGS1ClampSub1;
  G4RotationMatrix RotLGS1ClampSub2;
  RotLGS1ClampSub1.rotateY(0*deg);
  RotLGS1ClampSub2.rotateY(90*deg);

  solidLGS1Clamp = new G4Box("solidLGS1Clamp",0.5*50.0*mm,0.5*45.0*mm,0.5*26.0*mm);
  solidLGS1ClampSub3 = new G4Tubs("solidLGS1ClampSub1",0.5*0.0*mm,0.5*15.4*mm,26.2*mm,0*deg,360*deg);
  solidLGS1ClampSub2 = new G4Tubs("solidLGS1ClampSub2",0.5*0.0*mm,0.5*45.0*mm,10.0*mm,0*deg,360*deg);
  solidLGS1ClampSub1 = new G4Tubs("solidLGS1ClampSub3",120*mm,130.0*mm,51.0*mm,0*deg,360*deg);
  
  solidLGS1ClampFinal1 = new G4SubtractionSolid("solidLGS1ClampFinal1",solidLGS1Clamp,solidLGS1ClampSub1,
  						G4Transform3D(RotLGS1ClampSub2,G4ThreeVector(0,0,143.0*mm)));
  solidLGS1ClampFinal2 = new G4SubtractionSolid("solidLGS1ClampFinal2",solidLGS1ClampFinal1,solidLGS1ClampSub2,
  						G4Transform3D(RotLGS1ClampSub1,G4ThreeVector(0,0,-8.0*mm)));
  solidLGS1ClampFinal3 = new G4SubtractionSolid("solidLGS1ClampFinal2",solidLGS1ClampFinal2,solidLGS1ClampSub3,
  						G4Transform3D(RotLGS1ClampSub1,G4ThreeVector(0,0,0)));
  
  logicalLGS1ClampFinal = new G4LogicalVolume(solidLGS1ClampFinal3,Al,"solidLGS1ClampLV",0,0,0);
  

  G4RotationMatrix RotSolidClamp;
  
  G4double angleClamp=0;
  
  n=0;
  for(G4int j=0;j<NumberOfS1;j++)
    {
      G4double posClampx=(hightS1/cm+0.1)*(j-2)*cm;
      for(G4int i=0;i<3;i++)
	{
	  G4double posClampy=(130.0+0.5*26.0)*sin(angleClamp)*mm;
	  G4double posClampz=(-130.0-0.5*26.0)*cos(angleClamp)*mm;
	  
	  physiLGS1ClampFinal = new G4PVPlacement(G4Transform3D(RotSolidClamp,G4ThreeVector(posClampx,posClampy,posClampz)),
						  "solidLGS1ClampPV",logicalLGS1ClampFinal,physiWorld,false,n);
	  n++;
	  angleClamp += 120*deg;
	  RotSolidClamp.rotateX(120*deg);
	}
      G4double angleClamp2 = 360*deg/(NumberOfS1*3);
      angleClamp += angleClamp2;
      RotSolidClamp.rotateX(angleClamp2);
    }

  solidPMS1SH = new G4Tubs ("solidPMS1SH",0.5*38.1*mm,0.5*44.5*mm,0.5*160.0*mm,0*deg,360*deg);
  logicalPMS1SH = new G4LogicalVolume(solidPMS1SH,PM,"solidPMS1SHLV",0,0,0);
  
  G4double anglePMSH=0;
  G4RotationMatrix RotPMSH;
  RotPMSH.rotateX(180.0*deg);
  
  n=0;
  for(G4int j=0;j<NumberOfS1;j++)
    {
      G4double posPMSHx=(hightS1/cm+0.1)*(j-2)*cm;
      for(G4int i=0;i<3;i++)
	{
	  G4double posPMSHy=(144+0.5*160)*sin(anglePMSH)*mm;
	  G4double posPMSHz=(-144-0.5*160)*cos(anglePMSH)*mm;
	  physiPMS1SH = new G4PVPlacement(G4Transform3D(RotPMSH,G4ThreeVector(posPMSHx,posPMSHy,posPMSHz)),
					  "PMS1SHPV",logicalPMS1SH,physiWorld,false,n);
	  anglePMSH += 120*deg;
	  RotPMSH.rotateX(120*deg);
	}
      G4double anglePM2SH = 360*deg/(NumberOfS1*3);
      anglePMSH += anglePM2SH;
      RotPMSH.rotateX(anglePM2SH);
    }

  n=0;
  G4RotationMatrix RotPlate;
  RotPlate.rotateY(0.0*deg);
  G4RotationMatrix RotPlateHole;
  RotPlateHole.rotateY(90.0*deg);
  
  solidPlate = new G4Box("solidPlate",0.5*25*mm,0.5*1873.2*mm,0.5*1873.2*mm);
  //solidPlate = new G4Trd("solidPlate",1442.2*mm, 1442.2*mm, 1231.7*mm, 0.0*mm,0.5*25*mm);
  solidPlateHole = new G4Tubs("solidPlateHole",0.0*mm,0.5*150.0*mm,0.5*25.5*mm,0*deg,360*deg);
  solidPlateFinal = new G4SubtractionSolid("solidPlateFinal",solidPlate,solidPlateHole,
					   G4Transform3D(RotPlateHole,G4ThreeVector(0,0,0)));
  logicalPlate = new G4LogicalVolume(solidPlateFinal,Al,"solidPlateFinalLV",0,0,0);
  //physiPlate = new G4PVPlacement(G4Transform3D(RotPlate,G4ThreeVector(-62.5*mm,0.0*mm,190.5*mm)),
  //			 "solidPlatePV",logicalPlate,physiWorld,false,n);
  physiPlate = new G4PVPlacement(G4Transform3D(RotPlate,G4ThreeVector(-62.5*mm,0.0*mm,0.0*mm)),
  			 "solidPlatePV",logicalPlate,physiWorld,false,n);
  
  n=0;
  G4RotationMatrix RotS2Ring;
  RotS2Ring.rotateY(90.0*deg);
  
  solidS2Ring = new G4Tubs("solidS2Ring",700.0*mm,861.5*mm,0.5*25.0*mm,0*deg,360*deg);
  logicalS2Ring = new G4LogicalVolume(solidS2Ring,Al,"solidS2RingLV",0,0,0);
  physiS2Ring = new G4PVPlacement(G4Transform3D(RotS2Ring,G4ThreeVector(-50.0*mm+548.88*mm+0.5*25.0*mm,0,0)),
				  "solidS2RingPV",logicalS2Ring,physiWorld,false,n);
  
  n=0;
  G4RotationMatrix RotS1S2Leg;
  RotS1S2Leg.rotateX(0*deg);
  solidS1S2Leg = new G4Box("solidS1S2Leg",0.5*548.88*mm,0.5*100.0*mm,0.5*100.0*mm);
  solidS1S2LegSub = new G4Box("solidS1S2Leg",0.5*549.0*mm,0.5*94.0*mm,0.5*94.0*mm);
  solidS1S2LegFinal = new G4SubtractionSolid("solidS1S2LegFinal",solidS1S2Leg,solidS1S2LegSub,
					     G4Transform3D(RotS1S2Leg,G4ThreeVector(0,0,0)));
  logicalS1S2Leg = new G4LogicalVolume(solidS1S2LegFinal,Al,"solidS1S2LegFinalLV",0,0,0);
  G4double angleLeg = 0.0*deg;
  G4double angleLeg2 = 0.0*deg;
  for(G4int i=0;i<3;i++)
    {
      RotS1S2Leg.rotateX(angleLeg);
      G4double yPosLeg = 800.0*cos(angleLeg2);
      G4double zPosLeg = -800.0*sin(angleLeg2);
      physiS1S2Leg = new G4PVPlacement(G4Transform3D(RotS1S2Leg,G4ThreeVector(224.44*mm,yPosLeg,zPosLeg)),
				       "solidS1S2LegFinalPV",logicalS1S2Leg,physiWorld,false,n);
      angleLeg = 60.0*deg;
      angleLeg2 += 120.0*deg;
      n++;
    }
  
  n=0;
  G4RotationMatrix RotS2FootSub;
  RotS2FootSub.rotateX(0*deg);
  G4RotationMatrix RotS2Foot;
  RotS2Foot.rotateX(0*deg);
  G4RotationMatrix RotS2Plate;
  RotS2Plate.rotateZ(-35.0*deg);
  solidS2Foot = new G4Trap("solidS2Foot",100.0*mm,100.0*mm,197.0*mm,54.2*mm);
  solidS2FootSub = new G4Box("solidS2FootSub",0.5*197.1*mm,0.5*94.0*mm,0.5*94.0*mm);
  solidS2FootFinal = new G4SubtractionSolid("solidS2FootFinal",solidS2Foot,solidS2FootSub,
					    G4Transform3D(RotS2FootSub,G4ThreeVector(0,0,0)));
  logicalS2Foot = new G4LogicalVolume(solidS2FootFinal,Al,"solidS2FootLV",0,0,0);
  
  solidS2Plate = new G4Box("solidS2Plate",0.5*312.0*mm,0.5*8.0*mm,0.5*120.0*mm);
  logicalS2Plate = new G4LogicalVolume(solidS2Plate,Al,"solidS2PlateLV",0,0,0);
  
  solidS2Plate2 = new G4Box("solidS2Plate2",0.5*312.0*mm,0.5*8.0*mm,0.5*120.0*mm);
  logicalS2Plate2 = new G4LogicalVolume(solidS2Plate,Al,"solidS2Plate2LV",0,0,0);
  
  for(G4int j=0;j<NumberOfS2;j++)
    {
      G4double footangle = 360*deg/NumberOfS2*j;
      G4double yPosFoot = (805.0)*cos(footangle)*mm;
      G4double zPosFoot = -805.0*sin(footangle)*mm;
      G4double yPosPlate = (805.0-14.71)*cos(footangle)*mm;
      G4double zPosPlate = -805.0*sin(footangle)*mm;
      G4double yPosPlate2 = (805.0-7.0)*cos(footangle)*mm;
      G4double zPosPlate2 = -805.0*sin(footangle)*mm;
      physiS2Foot = new G4PVPlacement(G4Transform3D(RotS2Foot,G4ThreeVector(-62.5*mm+548.88*mm+0.5*197.0*mm,
									    yPosFoot,zPosFoot)),
				      "solidS2FootPV",logicalS2Foot,physiWorld,false,n);
      physiS2Plate = new G4PVPlacement(G4Transform3D(RotS2Plate,G4ThreeVector(-50.0*mm+548.88*mm+146.6*mm+25.0*mm+4.0*mm,
									      yPosPlate,zPosPlate)),
				       "solidS2PlatePV",logicalS2Plate,physiWorld,false,n);
      physiS2Plate2 = new G4PVPlacement(G4Transform3D(RotS2Plate,G4ThreeVector(-50.0*mm+548.88*mm+146.6*mm+25.0*mm+8.0*mm+2.0*mm+4.0*mm,
									       yPosPlate2,zPosPlate2)),
					"solidS2Plate2PV",logicalS2Plate2,physiWorld,false,n);
      G4double footangle2 = 360*deg/NumberOfS2;
      RotS2Foot.rotateX(-footangle2);
      RotS2Plate.rotateX(-footangle2);
      n++;
    }
  
  n=0;
  G4RotationMatrix RotS1Lid;
  RotS1Lid.rotateY(90.0*deg);
  solidS1Lid = new G4Tubs("solidS1Lid",75.0*mm,400*mm,0.5*4.0*mm,0*deg,360*deg);
  logicalS1Lid = new G4LogicalVolume(solidS1Lid,Al,"solidS1LidLV",0,0,0);
  physiS1Lid = new G4PVPlacement(G4Transform3D(RotS1Lid,G4ThreeVector(52*mm,0*mm,0*mm)),
				 "solidS1LidPV",logicalS1Lid,physiWorld,false,n);
  
  G4RotationMatrix RotS1LidS;
  RotS1LidS.rotateY(90.0*deg);
  solidS1LidS = new G4Tubs("solidS1LidS",397.0*mm,400*mm,0.5*100.0*mm,0*deg,360*deg);
  logicalS1LidS = new G4LogicalVolume(solidS1LidS,Al,"solidS1LidSLV",0,0,0);
  physiS1LidS = new G4PVPlacement(G4Transform3D(RotS1LidS,G4ThreeVector(0*mm,0*mm,0*mm)),
				 "solidS1LidSPV",logicalS1LidS,physiWorld,false,n);
  
  G4RotationMatrix RotS1Window;
  RotS1Window.rotateY(90.0*deg);
  solidS1Window = new G4Tubs("solidS1Window",0.0*mm,97.0*mm,0.5*0.1*mm,0*deg,360*deg);
  logicalS1Window = new G4LogicalVolume(solidS1Window,Al,"solidS1WindowLV",0,0,0);
  physiS1Window = new G4PVPlacement(G4Transform3D(RotS1Window,G4ThreeVector(54.0*mm+0.1*0.5*mm,0,0)),
				    "solidS1WindowPV",logicalS1Window,physiWorld,false,n);
  n++;
  physiS1Window = new G4PVPlacement(G4Transform3D(RotS1Window,G4ThreeVector(-75.0*mm+0.1*0.5*mm,0,0)),
				    "solidS1WindowPV",logicalS1Window,physiWorld,false,n);
  
  n=0;
  /**
     const G4int NUM = 14;
     G4double PPhoton[NUM] = {2.53*eV,2.58*eV,2.64*eV,2.70*eV,2.76*eV,2.82*eV,2.88*eV,
     2.95*eV,3.02*eV,3.10*eV,3.18*eV,3.26*eV,3.35*eV,3.44*eV};
     G4double RindexBC404[NUM] = {1.58,1.58,1.58,1.58,1.58,1.58,1.58,
     1.58,1.58,1.58,1.58,1.58,1.58,1.58};
     G4double RindexLG[NUM] = {1.49,1.49,1.49,1.49,1.49,1.49,1.49,
     1.49,1.49,1.49,1.49,1.49,1.49,1.49};
     G4double RindexAir[NUM] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0,1.0};
     
     G4double AbsBC404[NUM] = {140.0*cm,140.0*cm,140.0*cm,140.0*cm,140.0*cm,140.0*cm,140.0*cm,
     140.0*cm,140.0*cm,140.0*cm,140.0*cm,140.0*cm,140.0*cm,140.0*cm};
     G4double AbsLG[NUM] = {49.5*cm,49.5*cm,49.5*cm,49.5*cm,49.5*cm,49.5*cm,49.5*cm,
     49.5*cm,49.5*cm,49.5*cm,49.5*cm,49.5*cm,49.5*cm,49.5*cm};
     G4double scintillation[NUM] = {0.01,0.02,0.03,0.05,0.08,0.13,0.21,
     0.32,0.50,0.68,0.95,0.92,0.52,0.10};
     
     G4MaterialPropertiesTable *mptBC404 = new G4MaterialPropertiesTable();
     mptBC404 -> AddProperty("RINDEX",PPhoton,RindexBC404,NUM);
     mptBC404 -> AddProperty("ABSLENGTH",PPhoton,AbsBC404,NUM);
     mptBC404 -> AddProperty("SCINTILLATION",PPhoton,scintillation,NUM);
     BC404 -> SetMaterialPropertiesTable(mptBC404);
     
     G4MaterialPropertiesTable *mptLG = new G4MaterialPropertiesTable();
     mptLG -> AddProperty("RINDEX",PPhoton,RindexLG,NUM);
     mptLG -> AddProperty("ABSLENGTH",PPhoton,AbsLG,NUM);
     LG -> SetMaterialPropertiesTable(mptLG);
     
     
     G4MaterialPropertiesTable *mptAir = new G4MaterialPropertiesTable();
     mptAir -> AddProperty("RINDEX",PPhoton,RindexAir,NUM);
     Air -> SetMaterialPropertiesTable(mptAir);
     
     
     G4OpticalSurface *OpBC404Surface = new G4OpticalSurface("BC404Surface");
     //G4LogicalBorderSurface *BC404Surface = new G4LogicalBorderSurface("BC404Surface",physiS2,physiWorld,OpBC404Surface);
     //G4LogicalSkinSurface *BC404Surface = new G4LogicalSkinSurface("BC404Surface",logicalWorld,OpBC404Surface);
     OpBC404Surface->SetModel(unified);
     OpBC404Surface->SetType(dielectric_dielectric);
     OpBC404Surface->SetFinish(polished);
     
     //if(BC404Surface->GetVolume1()==physiS2) G4cout<<"Equal"<<G4endl;
     //if(BC404Surface->GetVolume2()==physiWorld) G4cout<<"Equal"<<G4endl;
     
     const G4int ENTRIS = 2;
     G4double PPhoton2[ENTRIS] = {2.53*eV,3.44*eV};
     G4double RindexBC404Air[ENTRIS] = {1.58,1.00};
     G4double specularlobeconst[ENTRIS] = {0.3,0.3};
     G4double specularspikeconst[ENTRIS] = {0.2,0.2};
     G4double backscatterconst[ENTRIS] = {0.1,0.1};
     


     G4MaterialPropertiesTable *BC404Air = new G4MaterialPropertiesTable();
     BC404Air -> AddProperty("RINDEX",PPhoton2,RindexBC404Air,ENTRIS);
     BC404Air -> AddProperty("SPECULARLOBECONSTANT",PPhoton2,specularlobeconst,ENTRIS);
     BC404Air -> AddProperty("SPECULARSPIKECONSTANT",PPhoton2,specularspikeconst,ENTRIS);
     BC404Air -> AddProperty("BACKSCATTERCONSTANT",PPhoton2,backscatterconst,ENTRIS);
     
     OpBC404Surface -> SetMaterialPropertiesTable(BC404Air);
     
     G4OpticalSurface *OpBC404LGS2Surface = new G4OpticalSurface("BC404LGS2Surface");
     //G4LogicalBorderSurface *BC404LGS2Surface = new G4LogicalBorderSurface("BC404LGS2Surface",physiS2,physiLGS2,OpBC404LGS2Surface);
     //G4LogicalSkinSurface *BC404LGS2Surface = new G4LogicalSkinSurface("BC404LGS2Surface",logicalLGS2,OpBC404LGS2Surface);
     OpBC404LGS2Surface->SetModel(unified);
     OpBC404LGS2Surface->SetType(dielectric_dielectric);
     OpBC404LGS2Surface->SetFinish(polished);
     
     G4double RindexBC404LGS2[ENTRIS] = {1.58,1.49};
     G4MaterialPropertiesTable *BC404LGS2 = new G4MaterialPropertiesTable();
     BC404LGS2 -> AddProperty("RINDEX",PPhoton2,RindexBC404LGS2,ENTRIS);
     BC404LGS2 -> AddProperty("SPECULARLOBECONSTANT",PPhoton2,specularlobeconst,ENTRIS);
     BC404LGS2 -> AddProperty("SPECULARSPIKECONSTANT",PPhoton2,specularspikeconst,ENTRIS);
     BC404LGS2 -> AddProperty("BACKSCATTERCONSTANT",PPhoton2,backscatterconst,ENTRIS);
     
     OpBC404LGS2Surface -> SetMaterialPropertiesTable(BC404LGS2);
  **/
  /**
  G4OpticalSurface *OpLGS2AirSurface = new G4OpticalSurface("LGS2AirSurface");
  G4LogicalBorderSurface *LGS2AirSurface = new G4LogicalBorderSurface("LGS2AirSurface",physiLGS2,physiWorld,OpLGS2AirSurface);
  OpLGS2AirSurface->SetModel(glisur);
  OpLGS2AirSurface->SetType(dielectric_dielectric);
  OpLGS2AirSurface->SetFinish(polished);
  
  G4double RindexLGS2Air[ENTRIS] = {1.49,1.00};
  G4MaterialPropertiesTable *LGS2Air = new G4MaterialPropertiesTable();
  LGS2Air -> AddProperty("RINDEX",PPhoton2,RindexLGS2Air,ENTRIS);
  LGS2Air -> AddProperty("SPECULARLOBECONSTANT",PPhoton2,specularlobeconst,ENTRIS);
  LGS2Air -> AddProperty("SPECULARSPIKECONSTANT",PPhoton2,specularspikeconst,ENTRIS);
  LGS2Air -> AddProperty("BACKSCATTERCONSTANT",PPhoton2,backscatterconst,ENTRIS);
  
  OpLGS2AirSurface-> SetMaterialPropertiesTable(LGS2Air);
  **/
  
  //---- Visualization attributes --------------------------------
  
  logicalWorld -> SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleTubeVisAtt = new G4VisAttributes(G4Colour(1.0,0.06,0.06));
  //G4VisAttributes* simpleTubeVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  simpleTubeVisAtt -> SetVisibility(true);
  logicalS1 -> SetVisAttributes(simpleTubeVisAtt);
  
  G4VisAttributes* simpleTrdVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  simpleTrdVisAtt -> SetVisibility(true);
  logicalS2 -> SetVisAttributes(simpleTrdVisAtt);

  G4VisAttributes* allLG = new G4VisAttributes(G4Colour(1.0,1.0,0.9));
  allLG -> SetVisibility(true);
  logicalLGS2  -> SetVisAttributes(allLG);
  logicalLGS2_b  -> SetVisAttributes(allLG);
  
  G4VisAttributes* allPM = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  allPM -> SetVisibility(true);
  logicalPMS1 -> SetVisAttributes(allPM);
  logicalPMS2 -> SetVisAttributes(allPM);
  
  G4VisAttributes *allCon = new G4VisAttributes(G4Colour(0.2,0.3,0.5));
  allCon -> SetVisibility(true);
  logicalLGS1Ring -> SetVisAttributes(allCon);
  
  G4VisAttributes *allClamp = new G4VisAttributes(G4Colour(0.5,0.3,0.5));
  allClamp -> SetVisibility(true);
  logicalLGS1ClampFinal -> SetVisAttributes(allClamp);
  
  G4VisAttributes *allPMSH = new G4VisAttributes(G4Colour(0.5,0.8,0.5));
  allPMSH -> SetVisibility(true);
  logicalPMS1SH -> SetVisAttributes(allPMSH);
  logicalPMS2SH -> SetVisAttributes(allPMSH);
  
  G4VisAttributes *allPlate = new G4VisAttributes(G4Colour(0.6,0.7,0.5));
  allPlate -> SetVisibility(true);
  logicalPlate -> SetVisAttributes(allPlate);
  
  G4VisAttributes *allS1S2Leg = new G4VisAttributes(G4Colour(0.7,0.6,0.5));
  allS1S2Leg -> SetVisibility(true);
  logicalS1S2Leg -> SetVisAttributes(allS1S2Leg);
  
  G4VisAttributes *allS2Ring = new G4VisAttributes(G4Colour(0.2,0.8,0.5));
  allS2Ring -> SetVisibility(true);
  logicalS2Ring -> SetVisAttributes(allS2Ring);
  
  G4VisAttributes *allS1Box = new G4VisAttributes(G4Colour(0.2,0.1,0.5));
  allS1Box -> SetVisibility(false);
  logicalS1Lid -> SetVisAttributes(allS1Box);
  logicalS1LidS -> SetVisAttributes(allS1Box);
  logicalS1Window -> SetVisAttributes(allS1Box);
  
  return physiWorld;
  
}

void S1S2DetCon::DefineMaterials()
{
  //---- Materials -------------------------------------
  G4double a;
  G4double z, iz;
  G4double density;
  G4String name,symbol;
  G4int nel, natoms;
  
  //----Element Carbon ---------------------------------
  a = 12.011*g/mole;
  iz = 6.;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", iz, a);
  //----Element Hydrogen -------------------------------
  a = 1.008*g/mole;
  iz = 1.;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz, a);
  //---- Composing Detector Material -------------------
  density = 1.032*g/cm3;
  BC404 = new G4Material(name="BC404", density, nel=2);
  BC404->AddElement(elC, natoms = 9);
  BC404->AddElement(elH, natoms = 10);
  
  //---- Element Nitrogrn ------------------------------
  a = 14.01*g/mole;
  iz = 7.;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz, a);
  //---- Element Oxigen --------------------------------
  a = 16.00*g/mole;
  iz = 8.;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz, a);
  //---- Element Iron ----------------------------------
  a = 55.85*g/mole;
  iz = 26.;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", iz, a);
  //---- Element Aluminium ----------------------------------
  a=26.98*g/mole;
  iz=13.;
  G4Element *elAl = new G4Element(name="Aluminium",symbol="Al", iz, a);

  //---- Composing Air ---------------------------------
  density = 1.29*mg/cm3;
  Air = new G4Material(name="Air", density, nel=2, kStateGas);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);
  
  density = 1.19*g/cm3;
  LG = new G4Material(name="LG",density,nel=3);
  LG->AddElement(elC, natoms = 5);
  LG->AddElement(elH, natoms = 8);
  LG->AddElement(elO, natoms = 2);
  
  //-- Composing PM-tubs of Iron -----------------------
  density = 7.87*g/cm3;
  PM = new G4Material(name="PM",density,nel=1);
  PM->AddElement(elFe,1.0);
  
  density = 2.70*g/cm3;
  Al = new G4Material(name="Al",density,nel=1);
  Al->AddElement(elAl,1.0);
			     
} 


G4double S1S2DetCon::GetxWorld()
{
    return xWorld;
}

G4double S1S2DetCon::GetoutRadS1()
{
    return outRadS1;
}

int S1S2DetCon::GetNumberOfS1()
{
    return NumberOfS1;
}

int S1S2DetCon::GetNumberOfS2()
{
    return NumberOfS2;
}

void S1S2DetCon::SetxTrd(G4double value)
{
  x1Trd = value;
  x2Trd = value;
}

void S1S2DetCon::Sety1Trd(G4double value)
{
  y1Trd =  value;
}

void S1S2DetCon::Sety2Trd(G4double value)
{
  y2Trd =  value;
}

void S1S2DetCon::SetzTrd(G4double value)
{
  zTrd =  value;
}

void S1S2DetCon::SetxposS2Abs(G4double value)
{
  xposS2Abs =  value;
}

void S1S2DetCon::SetyposS2Abs(G4double value)
{
  yposS2Abs =  value;
}

void S1S2DetCon::SetzposS2Abs(G4double value)
{
  zposS2Abs =  value;
}

void S1S2DetCon::Setangle1(G4double value)
{
  meanSangle = value;
  angle1=2*meanSangle-tiltangle;
  
  rangle = (180*deg-2*meanSangle);
  meanFlightpath = sqrt(2*pow(704.6,2)*(1-cos(rangle)))*mm;
  
  xposS2Abs = meanFlightpath*cos(meanSangle);
  yposS2Abs = meanFlightpath*sin(meanSangle);
  zposS2Abs = meanFlightpath*sin(meanSangle);
  
  //y1Trd = 2*(yposS2Abs-zTrd/2*cos(angle1))*3.1416/NumberOfS2-0.5*cm;
  //y2Trd = 2*(yposS2Abs+zTrd/2*cos(angle1))*3.1416/NumberOfS2-0.5*cm;
  y1Trd = 95.0*mm;
  y2Trd = 134.7*mm;
  G4cout<<"y1Trd: "<<y1Trd<<G4endl;
  G4cout<<"y2Trd: "<<y2Trd<<G4endl;
}

G4double S1S2DetCon::Getangle1()
{
  return angle1;
}

G4double S1S2DetCon::GetFlightpath()
{
  return meanFlightpath;
}

G4double S1S2DetCon::GetMeanSangle()
{
  return meanSangle;
}

void S1S2DetCon::SetS1Thick(G4double value)
{
  hightS1 = value;
}

void S1S2DetCon::SetNumberOfS2(G4int value)
{
  NumberOfS2 = value;
}

void S1S2DetCon::SetNumberOfS1(G4int value)
{
  NumberOfS1 = value;
}

void S1S2DetCon::SetS1Radius(G4double value)
{
  outRadS1=value;
}

void S1S2DetCon::SetTiltangle(G4double value)
{
  tiltangle = value;
}

void S1S2DetCon::UpdateGeometry()
{
  G4RunManager::GetRunManager() -> DefineWorldVolume(ConstructAll());
}
