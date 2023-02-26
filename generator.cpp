#include "generator.h"
#include "gmsh.h"

#include "element2d.h"
#include "czinsertiontool2d.h"
#include "cohesivezone2d.h"

#include <spdlog/spdlog.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_set>



void Generator::LoadFromFileWithCrop(std::string MSHFileName, double indenterX, double indenterY, int cropType)
{
    spdlog::info("loading from file {}",MSHFileName);
    gmsh::initialize();
    gmsh::open(MSHFileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, int> mtags; // gmsh nodeTag -> node index

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue;
        icy::Node2D *nd = mesh2d.AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = Eigen::Vector2d(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // GET ELEMENTS - per grain (entity)
    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,2);

    std::unordered_set<int> used_nodes;

    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> trisTags, nodeTagsInTris;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris,entityTag);

        for(std::size_t i=0;i<trisTags.size();i++)
        {
            icy::Node2D* nds[3];
            bool crop = false;
            for(int k=0;k<3;k++)
            {
                nds[k] = &mesh2d.nodes[mtags.at(nodeTagsInTris[i*3+k])];
                Eigen::Vector2d &ndpos = nds[k]->x0;
                double x = ndpos.x();
                double y = ndpos.y();

                // skip the elements, whcih need to be cropped
                // cropType 1 -> normal block with indenter at 1.0
                // cropType 2 -> normal block with indenter at 2.0
                // cropType 3 -> tall block with indenter at 1.0
                // cropType 4 -> tall block with indenter at 2.0
                // cropType 5 -> small keel with indenter at 5.0
                // cropType 6 -> small keel with indenter at 6.0
                // cropType 7 -> tall keel with indenter at 5.0
                // cropType 8 -> tall keel with indenter at 6.0
                if(cropType >=1 && cropType <=4 && x > indenterX && y > (indenterY-indenterRadius)) crop = true;


            }
            if(crop) continue;

            icy::Element2D *elem = mesh2d.AddElement();
            elem->grainId = (int)j;
            for(int k=0;k<3;k++)
            {
                elem->nds[k] = nds[k]->globId;
                used_nodes.insert(nds[k]->globId);
            }
        }
    }
    mesh2d.nodes.erase(std::remove_if(mesh2d.nodes.begin(), mesh2d.nodes.end(),
                                      [used_nodes](const icy::Node2D &nd){ return used_nodes.count(nd.globId)==0 ? true : false;}),
            mesh2d.nodes.end());

    std::unordered_map<int,int> updatedNodeIds; // oldId -> newId

    for(int i=0;i<mesh2d.nodes.size();i++)
    {
        updatedNodeIds.insert({mesh2d.nodes[i].globId,i});
        mesh2d.nodes[i].globId = i;
    }

    for(icy::Element2D &elem : mesh2d.elems)
        for(int k=0;k<3;k++) elem.nds[k] = updatedNodeIds.at(elem.nds[k]);

    gmsh::finalize();

    icy::CZInsertionTool2D czit;
    czit.InsertCZs(mesh2d);

    for(icy::Element2D &elem : mesh2d.elems) elem.Precompute(mesh2d.nodes);     // Dm matrix and volume


    // attach bottom
    if(cropType >=1 && cropType <=4)
    {
        for(icy::Node2D &nd : mesh2d.nodes)
            if(nd.x0.y()==0) nd.group = 2;
    }
    else if(cropType == 5 || cropType == 6)
    {
        for(icy::Node2D &nd : mesh2d.nodes)
            if(nd.x0.y()>1.0-1e-5) nd.group = 2;
    }
    else
    {
        for(icy::Node2D &nd : mesh2d.nodes)
            if(nd.x0.y()>2.0-1e-5) nd.group = 2;
    }


    spdlog::info("nds {}; elems {}; czs {}; grains {}", mesh2d.nodes.size(), mesh2d.elems.size(), mesh2d.czs.size(), dimTagsGrains.size());

    if(!QDir("results").exists())
        QDir().mkdir("results");

    QFileInfo fi(QString(MSHFileName.c_str()));

    CreatePy2D(fi.baseName().toStdString());

    spdlog::info("LoadFromFile done");
}


void Generator::CreatePy2D(std::string outputFileName)
{
    // SAVE as .py
    std::ofstream s;
    s.open(outputFileName+".py", std::ios_base::trunc|std::ios_base::out);
    s << std::setprecision(9);
    s << "from abaqus import *\n";
    s << "from abaqusConstants import *\n";
    s << "from caeModules import *\n";

    s << "import mesh\n";
    s << "import regionToolset\n";
    s << "import os\n";

    s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)\n";

    for(icy::Node2D *nd : mesh2d.nodes)
        s << "p.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << ",0))\n";

    s << "n = p.nodes\n";

    for(icy::Element2D *e : mesh2d.elems)
        s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<e->nds[1]->globId<<
             "],n["<<e->nds[2]->globId<<"]), elemShape=TRI3)\n";

    for(icy::CohesiveZone2D *c : mesh2d.czs)
        s << "p.Element(nodes=(n["<<c->nds[0]->globId<<
             "], n["<<c->nds[1]->globId<<
             "], n["<<c->nds[3]->globId<<
             "], n["<<c->nds[2]->globId<<
             "]), elemShape=QUAD4)\n";

    s << "elemType_bulk = mesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD, secondOrderAccuracy=OFF,"
         "distortionControl=ON, lengthRatio=0.1, elemDeletion=ON)\n";

    bool hasCZs = mesh2d.czs.size()>0;

    if(hasCZs)
        s << "elemType_coh = mesh.ElemType(elemCode=COH2D4, elemLibrary=STANDARD,"
             "distortionControl=ON, lengthRatio=0.1, elemDeletion=ON)\n";


    // region1 - bulk elements
    s << "region1 = p.elements[0:" << mesh2d.elems.size() << "]\n";
    s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
    s << "p.Set(elements=(region1,), name='Set-1-elems')\n";

    if(hasCZs)
    {
        s << "region2cz = p.elements[" << mesh2d.elems.size() << ":" << mesh2d.elems.size() + mesh2d.czs.size() << "]\n";
        s << "p.setElementType(regions=(region2cz,), elemTypes=(elemType_coh,))\n";
        s << "p.Set(elements=(region2cz,), name='Set-2-czs')\n";
    }

    // region - pinned nodes

    s << "region3pinned = (";
    for(icy::Node2D *nd : mesh2d.nodes)
        if(nd->group==2)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";

    // region - all nodes
    s << "region4allNodes = (";
    for(icy::Node2D *nd : mesh2d.nodes)
    {
        if(nd->x0.y() < 0.8) continue; // only include the top layer
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    }
    s << ")\n";
    s << "p.Set(nodes=region4allNodes,name='Set4-all')\n";


    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((900.0, ), ))\n";
    s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

    CreateCDP(s);

    // plastic material
//        s << "mat1.Plastic(scaleStress=None,table="
//                "((50000.0, 0.0), (100000.0, 0.01), (1000000.0, 0.05), (10000000.0,0.1)))\n";

    // cz material
    if(hasCZs)
    {
        s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-czs')\n";
        s << "mat2.Density(table=((1.0, ), ))\n";
        s << "mat2.MaxsDamageInitiation(table=((" << czsStrength << "," << czsStrength*2 << "," << czsStrength*2 << "), ))\n";
        s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, softening=EXPONENTIAL, table=((" << czEnergy << ", ), ))\n";
        s << "mat2.Elastic(type=TRACTION, table=((" << czElasticity << "," << czElasticity << "," << czElasticity << "), ))\n";

        s << "mdb.models['Model-1'].CohesiveSection(name='Section-2-czs', "
             "material='Material-2-czs', response=TRACTION_SEPARATION, "
             "outOfPlaneThickness=None)\n";
    }

    // sections
    s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1-bulk', "
         "material='Material-1-bulk', thickness=None)\n";

    // section assignments
    s << "region = p.sets['Set-1-elems']\n";
    s << "p.SectionAssignment(region=region, sectionName='Section-1-bulk', offset=0.0, "
         "offsetType=MIDDLE_SURFACE, offsetField='', "
         "thicknessAssignment=FROM_SECTION)\n";

    if(hasCZs)
    {
        s << "region = p.sets['Set-2-czs']\n";
        s << "p = mdb.models['Model-1'].parts['MyPart1']\n";
        s << "p.SectionAssignment(region=region, sectionName='Section-2-czs', offset=0.0, "
             "offsetType=MIDDLE_SURFACE, offsetField='', "
             "thicknessAssignment=FROM_SECTION)\n";
    }

    // indenter
    s << "s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)\n";
    s << "g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n";
    s << "s.ArcByCenterEnds(center=(0.0, 0.0), point1=("<< -indenterRadius << ", 0.0), point2=(" << indenterRadius << ", -0.0125), direction=COUNTERCLOCKWISE)\n";

    s << "p2 = mdb.models['Model-1'].Part(name='Part-2', dimensionality=TWO_D_PLANAR, type=ANALYTIC_RIGID_SURFACE)\n";

    s << "p2.AnalyticRigidSurf2DPlanar(sketch=s)\n";


    s << "v1 = p2.vertices\n";

    s << "p2.ReferencePoint(point=p2.InterestingPoint(p2.edges[0], CENTER))\n";


    // assembly
    s << "a1 = mdb.models['Model-1'].rootAssembly\n";
    s << "a1.DatumCsysByDefault(CARTESIAN)\n";

    // add and rotate main part
    s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";

    const double &R = indenterRadius;
    const double &d = indentationDepth;
    double b = sqrt(R*R - (R-d)*(R-d));
    double xOffset = indenterOffset == 0 ? notchOffset-b : indenterOffset;
    xOffset -= interactionRadius*1.5;
    spdlog::info("xOffset {}; indenterOffset {}; notchOffset {}",xOffset, indenterOffset, notchOffset);
    //horizontalOffset += -sqrt(pow(indenterRadius,2)-pow(indenterRadius-indenterDepth,2))-1e-7;
    double yOffset = blockHeight-d+R;

    double zOffset = 0;
    s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
    // rotate indenter
    s << "a1.rotate(instanceList=('Part-2-1', ), axisPoint=(0.0, 0.0, 0.0)," << "axisDirection=(0.0, 0.0, 1.0), angle=45.0)\n";
    s << "a1.translate(instanceList=('Part-2-1', ), vector=("<< xOffset << ", " << yOffset << ", " << zOffset << "))\n";

    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timeToRun << ", improvedDtMethod=ON)\n";

    // create field output request
    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << nFrames <<
         ",variables=('S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', "
             "'U', 'V', 'A', 'RF', 'CSTRESS', 'DAMAGEC', 'DAMAGET', 'DAMAGESHR', 'EVF', "
             "'STATUS', 'SDEG'))\n";

    // gravity load
    s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

    // BC - pinned nodes
    s << "region = inst1.sets['Set3-pinned']\n";
    s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";

    // BC - moving indenter
    double xVelocity = indentationRate;
    double yVelocity = 0;
    s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
    s << "refPoints1=(r1[2], )\n";
    s << "region = a1.Set(referencePoints=refPoints1, name='Set-1-indenterRP')\n";
    s << "mdb.models['Model-1'].VelocityBC(name='BC-2', createStepName='Step-1', "
         "region=region, v1="<< xVelocity << ", v2=" << yVelocity << ", v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, "
                                                                     "amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')\n";

    s << "ed1 = a1.instances['Part-2-1'].edges\n";
    s << "side2Edges1 = ed1[0:1]\n";
    s << "a1.Surface(name='Surf-1', side2Edges=side2Edges1)\n";

    s << "mdb.models['Model-1'].RigidBody(name='Constraint-1', refPointRegion=regionToolset.Region(referencePoints="
         "(a1.instances['Part-2-1'].referencePoints[2],)), surfaceRegion=a1.surfaces['Surf-1'])\n";


    // rigid body constraint

    // create interaction property - hard collisions
    s << "mdb.models['Model-1'].ContactProperty('IntProp-2czs')\n";

    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].TangentialBehavior("
         "formulation=FRICTIONLESS)\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].NormalBehavior("
         "pressureOverclosure=HARD, allowSeparation=ON, "
         "constraintEnforcementMethod=DEFAULT)\n";

    // interaction property - penalty-based
    s << "mdb.models['Model-1'].ContactProperty('IntProp-1')\n";
      s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior("
        "formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, "
        "pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=(("
        "0.1, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, "
        "fraction=0.005, elasticSlipStiffness=None)\n";

       s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior("
        "pressureOverclosure=EXPONENTIAL, table=((20000000.0, 0.0), (0.0, " << interactionRadius << ")), "
        "maxStiffness=None, constraintEnforcementMethod=DEFAULT)\n";


    // create interaction itself
       if(createSelfCollisions)
       {
           s << "mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Step-1')\n";
           s << "mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep("
                "stepName='Step-1', useAllstar=ON)\n";
           s << "mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep("
                "stepName='Step-1', assignments=((GLOBAL, SELF, 'IntProp-2czs'), ))\n";
       }

    // additional contact between surface and nodes
    s << "mdb.models['Model-1'].SurfaceToSurfaceContactExp(clearanceRegion=None,"
        "createStepName='Initial', datumAxis=None, initialClearance=OMIT, "
        "interactionProperty='IntProp-1', main="
        "a1.surfaces['Surf-1'], "
        "mechanicalConstraint=KINEMATIC, name='Int-2', secondary="
        "a1.sets['MyPart1-1.Set4-all'], sliding=FINITE)\n";


    // record indenter force

    s << "mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name="
         "'H-Output-2', rebar=EXCLUDE, region="
         "mdb.models['Model-1'].rootAssembly.sets['Set-1-indenterRP'], sectionPoints="
         "DEFAULT, timeInterval=0.0001, variables=('RF1','RF2', ))\n";

    //create job
    s << "mdb.Job(name='" << outputFileName << "', model='Model-1', description='', type=ANALYSIS,"
                                        "atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,"
                                        "memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE,"
                                        "nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,"
                                        "contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',"
                                        "resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains="
      <<numberOfCores<<","
                       "activateLoadBalancing=False, numThreadsPerMpiProcess=1,"
                       "multiprocessingMode=DEFAULT, numCpus="<<numberOfCores<<")\n";

    s.close();
}

void Generator::CreateCDP(std::ofstream &s)
{
    if(!createCDP) return;

    std::vector<double> CompressionStrainTable
    {
        0,
        0.000241232889562297,
        0.000354471872727266,
        0.000521980525252549,
        0.00071138618855219,
        0.000933655856565658,
        0.00118044919730642,
        0.00149661554545456,
        0.00194945956969697,
        0.00256248560673401,
        0.00333461398720537,
        0.00438538037710438,
        0.00896666666188552
    };

    std::vector<double> ConcreteCompressionHardening {
        27887324,
        37126760,
        38535212,
        39718308,
        40000000,
        39887324,
        39267604,
        38028168,
        36000000,
        32901408,
        28957748,
        24225352,
        400000
    };

    std::vector<double> ConcreteCompressionDamage {
        0,
        0,
        0,
        0,
        0,
        0.00281690000000001,
        0.0183099,
        0.0492958,
        0.1,
        0.1774648,
        0.2760563,
        0.3943662,
        0.99
    };



    // TENSILE PARAMETERS
    std::vector<double> TensionStrainTable {
        0,
        4.82465779124578E-05,
        7.08943745454513E-05,
        0.000104396105050508,
        0.000142277237710436,
        0.00018673117131313,
        0.000236089839461282,
        0.00029932310909091,
        0.000389891913939393,
        0.000512497121346801,
        0.000666922797441073,
        0.000877076075420875,
        0.0017933333323771
    };

    std::vector<double> ConcreteTensionStiffening {
        5577464.8,
        7425352,
        7707042.4,
        7943661.6,
        8000000,
        7977464.8,
        7853520.8,
        7605633.6,
        7200000,
        6580281.6,
        5791549.6,
        4845070.4,
        80000
    };

    std::vector<double> ConcreteTensionDamage {
        0,
        0,
        0,
        0,
        0,
        0.00281690000000001,
        0.0183099,
        0.0492958,
        0.1,
        0.1774648,
        0.2760563,
        0.3943662,
        0.99
    };


    s << "mat1.ConcreteDamagedPlasticity(table=((40.0, 0.1, 1.16, 0.6667, 0.0), ))\n";

    s << std::setprecision(15);

    s << "mat1.concreteDamagedPlasticity.ConcreteCompressionHardening(table=(";
    for(int i=0;i<ConcreteCompressionHardening.size();i++)
    {
        const double stress = ConcreteCompressionHardening[i];
        const double strain = CompressionStrainTable[i];
        s << "(" << stress << ","  << strain << ")";
        if(i!=ConcreteCompressionHardening.size()-1) s << ",";
    }
    s << "))\n";

    s << "mat1.concreteDamagedPlasticity.ConcreteTensionStiffening(table=(";
    for(int i=0;i<ConcreteTensionStiffening.size();i++)
    {
        const double stress = ConcreteTensionStiffening[i];
        const double strain = TensionStrainTable[i];
        s << "(" << stress << ","  << strain << ")";
        if(i!=ConcreteTensionStiffening.size()-1) s << ",";
    }
    s << "))\n";

    s << "mat1.concreteDamagedPlasticity.ConcreteCompressionDamage(table=(";
    for(int i=0;i<ConcreteCompressionDamage.size();i++)
    {
        const double damage = ConcreteCompressionDamage[i];
        const double strain = CompressionStrainTable[i];
        s << "(" << damage << ","  << strain << ")";
        if(i!=ConcreteCompressionDamage.size()-1) s << ",";
    }
    s << "))\n";

    s << "mat1.concreteDamagedPlasticity.ConcreteTensionDamage(table=(";
    for(int i=0;i<ConcreteTensionDamage.size();i++)
    {
        const double damage = ConcreteTensionDamage[i];
        const double strain = TensionStrainTable[i];
        s << "(" << damage << ","  << strain << ")";
        if(i!=ConcreteTensionDamage.size()-1) s << ",";
    }
    s << "))\n";

    //*Concrete Failure,TYPE=Strain
    //0.0128889,  0.00348889, 0.9,  0.99

    std::cout << "*Concrete Failure,TYPE=Strain\n";
    std::cout << TensionStrainTable.back() << ",";
    std::cout << CompressionStrainTable.back() << ",";
    std::cout << ConcreteTensionDamage.back() << ",";
    std::cout << ConcreteCompressionDamage.back() << '\n';
}


