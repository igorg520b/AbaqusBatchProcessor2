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
    setupType = cropType;
    this->indenterY = indenterY;
    this->indenterX = indenterX;

    spdlog::info("loading from file {}",MSHFileName);
    gmsh::initialize();
    gmsh::open(MSHFileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, int> mtags; // gmsh nodeTag -> node index

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    spdlog::info("reading nodes");
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

    spdlog::info("reading per-grain elements");
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
                if(cropType >=1 && cropType <=4 && x < indenterX && y > (indenterY-indenterRadius)) crop = true;
                if(cropType >=5 && cropType <=8 && x > indenterX && y < (indenterRadius)) crop = true;


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
    const double iR = indenterRadius+1e-3;
    if(cropType >=1 && cropType <=4)
    {
        double alpha = forceAngle * pi/180.;
        fX = forceMagnitude * cos(alpha);
        fY = -forceMagnitude * sin(alpha);
        const double bH = (cropType == 1 || cropType == 2) ? 1.0 : 2.0; // block height
        double a1 = -1.570796326794897;
        double a2 = -1.16698;
        for(icy::Node2D &nd : mesh2d.nodes)
        {
            double y = nd.x0.y();
            double dx = nd.x0.x()-indenterX;
            double dy = nd.x0.y()-indenterY;
            if(y==0)
            {
                nd.group = 2;
            }
            else if(dx*dx+dy*dy <= iR*iR)
            {
                double a = atan2(dy,dx);
                if(a>=a1 && a<=a2) nd.group=5;
            }
        }

    }
    else if(cropType == 5 || cropType == 6)
    {
        double a1 = 1.570796326794897;
        double a2 = 2.356194490192345;
        double avg = (a1+a2)/2;
        fX = 1e6 * cos(avg);
        fY = 1e6 * sin(avg);
        for(icy::Node2D &nd : mesh2d.nodes)
        {
            double y = nd.x0.y();
            double dx = nd.x0.x()-indenterX;
            double dy = nd.x0.y()-indenterY;
            if(y>1.0-1e-5) nd.group = 2;
            else if(dx*dx+dy*dy <= iR*iR)
            {
                double a = atan2(dy,dx);
                if(a>=a1 && a<=a2) nd.group=5;
            }
        }
    }
    else if(cropType == 7 || cropType == 8)
    {
        double a1 = 1.570796326794897;
        double a2 = 2.356194490192345;
        double avg = (a1+a2)/2;
        fX = 1e6 * cos(avg);
        fY = 1e6 * sin(avg);
        for(icy::Node2D &nd : mesh2d.nodes)
        {
            double y = nd.x0.y();
            double dx = nd.x0.x()-indenterX;
            double dy = nd.x0.y()-indenterY;
            if(y>2.0-1e-5) nd.group = 2;
            else if(dx*dx+dy*dy <= iR*iR)
            {
                double a = atan2(dy,dx);
                if(a>=a1 && a<=a2) nd.group=5;
            }
        }
    }


    spdlog::info("nds {}; elems {}; czs {}; grains {}", mesh2d.nodes.size(), mesh2d.elems.size(), mesh2d.czs.size(), dimTagsGrains.size());



    QFileInfo fi(QString(MSHFileName.c_str()));
    fileName = fi.baseName().toStdString();

    CreatePy2D();
    CreateExportScript();

    spdlog::info("LoadFromFile done");
}


void Generator::CreatePy2D()
{
    // SAVE as .py
    std::ofstream s;
    s.open("results\\"+fileName+".py", std::ios_base::trunc|std::ios_base::out);
    s << std::setprecision(9);
    s << "from abaqus import *\n";
    s << "from abaqusConstants import *\n";
    s << "from caeModules import *\n";

    s << "import mesh\n";
    s << "import regionToolset\n";
    s << "import os\n";

    s << "os.chdir(r\"";
    s << QDir::currentPath().toStdString() << "\\inp";
    s << "\")\n";

    s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)\n";

    s << "print(\"importing " << mesh2d.nodes.size() << " nodes\")\n";
    for(const icy::Node2D &nd : mesh2d.nodes)
    {
        int k = nd.globId;
        if(k%(mesh2d.nodes.size()/100)==0)
            s << "print(\"" << 100*k/mesh2d.nodes.size() << "%\")\n";
        s << "p.Node(coordinates=(" << nd.x0[0] << "," << nd.x0[1] << ",0))\n";
    }

    s << "n = p.nodes\n";

    s << "print(\"importing " << mesh2d.elems.size() << " elements\")\n";
    for(const icy::Element2D &e : mesh2d.elems)
    {
        int k = e.elemId;
        if(k%(mesh2d.elems.size()/100)==0)
            s << "print(\"" << 100*k/mesh2d.elems.size() << "%\")\n";
        s << "p.Element(nodes=(n["<<e.nds[0]<<"],n["<<e.nds[1]<<"],n["<<e.nds[2]<<"]), elemShape=TRI3)\n";
    }

    s << "print(\"importing " << mesh2d.czs.size() << " czs\")\n";
    for(const icy::CohesiveZone2D &c : mesh2d.czs)
        s << "p.Element(nodes=(n["<<c.nds[0]<<"], n["<<c.nds[1]<<"], n["<<c.nds[3]<<"], n["<<c.nds[2]<<"]), elemShape=QUAD4)\n";

    s << "print(\"finalizing\")\n";

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
    for(const icy::Node2D &nd : mesh2d.nodes)
        if(nd.group==2)
            s << "p.nodes["<<nd.globId<<":"<<nd.globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";

    // forced nodes
    s << "region5forced = (";
    for(const icy::Node2D &nd : mesh2d.nodes)
        if(nd.group==5)
            s << "p.nodes["<<nd.globId<<":"<<nd.globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region5forced,name='Set5-forced')\n";

    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((900.0, ), ))\n";
    s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

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


    // assembly
    s << "a1 = mdb.models['Model-1'].rootAssembly\n";
    s << "a1.DatumCsysByDefault(CARTESIAN)\n";

    // add and rotate main part
    s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";

    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timeToRun << ", improvedDtMethod=ON)\n";

    // create field output request
    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << nFrames <<
         ",variables=('S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', "
             "'U', 'V', 'A', 'CSTRESS', 'DAMAGEC', 'DAMAGET', 'DAMAGESHR', 'EVF', "
             "'STATUS', 'SDEG'))\n";

    // gravity load
    s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

    // force load
    s << "mdb.models['Model-1'].TabularAmplitude(data=((0.0, 0.0), (0.1, 1.0)), name="
         "'Amp-1', smooth=SOLVER_DEFAULT, timeSpan=STEP)\n";

    s << "mdb.models['Model-1'].ConcentratedForce(amplitude='Amp-1', cf1=" << fX << ", cf2=" << fY <<
         ", createStepName='Step-1', distributionType=UNIFORM, field='', "
         "localCsys=None, name='Load-2', region="
         "mdb.models['Model-1'].rootAssembly.instances['MyPart1-1'].sets['Set5-forced'])\n";

    // BC - pinned nodes
    s << "region = inst1.sets['Set3-pinned']\n";
    s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";


    // rigid body constraint

    // create interaction property - hard collisions
    s << "mdb.models['Model-1'].ContactProperty('IntProp-2czs')\n";

    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].TangentialBehavior("
         "formulation=FRICTIONLESS)\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].NormalBehavior("
         "pressureOverclosure=HARD, allowSeparation=ON, "
         "constraintEnforcementMethod=DEFAULT)\n";


    //create job
    s << "mdb.Job(name='" << fileName << "', model='Model-1', description='', type=ANALYSIS,"
                                        "atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,"
                                        "memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE,"
                                        "nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,"
                                        "contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',"
                                        "resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains="
      <<numberOfCores<<","
                       "activateLoadBalancing=False, numThreadsPerMpiProcess=1,"
                       "multiprocessingMode=DEFAULT, numCpus="<<numberOfCores<<")\n";


     // write .inp file
    s << "mdb.jobs['" << fileName << "'].writeInput(consistencyChecking=OFF)\n";

    // write .cae
    s << "mdb.saveAs(\"" << QDir::currentPath().toStdString() + "/results/"+fileName+".cae\"" << ")\n";

    s.close();
}



void Generator::AddEntryToJobGeneratorBat(std::ofstream &s)
{
    s << AbaqusExecutablePath << " cae noGUI=\"";
    s << QDir::currentPath().toStdString() + "\\results\\"+fileName+".py\"\n";
}

void Generator::AddEntryToJobExecutionBat(std::ofstream &s)
{
    s << AbaqusExecutablePath << " job=" << fileName;
    s << " double=both output_precision=full cpus=12 interactive\n";
}

void Generator::AddEntryToBinaryExportBat(std::ofstream &s)
{
    s << AbaqusExecutablePath << " cae noGUI=\"";
    s << QDir::currentPath().toStdString() + "\\exportscripts\\"+fileName+".py\"\n";
}

void Generator::CreateExportScript()
{
    std::ofstream s;
    s.open("exportscripts\\"+fileName+".py", std::ios_base::trunc|std::ios_base::out);


    s << "from odbAccess import *\n"
    "from abaqusConstants import *\n"
    "from odbMaterial import *\n"
    "from odbSection import *\n"
    "from abaqus import *\n"
    "from abaqusConstants import *\n"
    "from caeModules import *\n"
    "import mesh\n"
    "import regionToolset\n"
    "import os\n"
    "import array\n"
    "import struct\n";

    s << "fileName='" << fileName << "'\n";
    s << "fodb = openOdb(path=r'" << QDir::currentPath().toStdString() << "/inp/'+fileName+'.odb')\n";
    s << "ma = fodb.rootAssembly\n"
    "f = open(fileName+'.bin', \"wb\")\n"
    "inst_block = ma.instances['MYPART1-1']\n"
    "nds = inst_block.nodes\n"
    "elems = inst_block.elements\n"
    "nNodes = len(nds)\n"
    "nElems = len(elems)\n\n"

    "nFrames = len(fodb.steps['Step-1'].frames)\n"
    "nNodeValues = len(fodb.steps['Step-1'].frames[0].fieldOutputs['U'].values)\n"
    "nStressValues = len(fodb.steps['Step-1'].frames[0].fieldOutputs['S'].values)\n"
    "nSDEGValues = len(fodb.steps['Step-1'].frames[0].fieldOutputs['SDEG'].values)\n"
    "nTrianglesInFrame = len(fodb.steps['Step-1'].frames[0].fieldOutputs['STATUS'].values)\n"

    "f.write(struct.pack('i', nNodes))\n"
    "f.write(struct.pack('i', nElems))\n"
    "f.write(struct.pack('i', nFrames))\n"
    "f.write(struct.pack('i', nNodeValues))\n"
    "f.write(struct.pack('i', nStressValues))\n"
    "f.write(struct.pack('i', nSDEGValues))\n"
    "f.write(struct.pack('i', nTrianglesInFrame))\n";
    s<< "f.write(struct.pack('i', " << setupType << "))\n";
    s<< std::setprecision(16);
    s<< "f.write(struct.pack('d', " << indenterX << "))\n";
    s<< "f.write(struct.pack('d', " << indenterY << "))\n";

    s << "arrNdLabels = array.array('i',(0 for i in range(0,nNodes)))\n"
    "arrNdCoordinates = array.array('d', (0 for i in range(0,nNodes*2)))\n"
    "for k in range(nNodes):\n"
    "    arrNdLabels[k] = nds[k].label\n"
    "    arrNdCoordinates[k*2+0] = nds[k].coordinates[0]\n"
    "    arrNdCoordinates[k*2+1] = nds[k].coordinates[1]\n"
    "buf = struct.pack('%si' % len(arrNdLabels), *arrNdLabels)\n"
    "f.write(buf)\n"
    "buf = struct.pack('%sd' % len(arrNdCoordinates), *arrNdCoordinates)\n"
    "f.write(buf)\n";


    s << "arrElemData = array.array('i', (0 for i in range(0,nElems*6)))\n"
    "nCZS = 0\n"
    "nTriangles = 0\n"
    "for k in range(nElems):\n"
    "    arrElemData[k*6+1] = elems[k].label\n"
    "    if elems[k].type == 'CPS3':\n"
    "        arrElemData[k*6+0] = 77 #type is triangle\n"
    "        arrElemData[k*6+2] = elems[k].connectivity[0]\n"
    "        arrElemData[k*6+3] = elems[k].connectivity[1]\n"
    "        arrElemData[k*6+4] = elems[k].connectivity[2]\n"
    "        arrElemData[k*6+5] = -1\n"
    "        nTriangles += 1\n"
    "    elif elems[k].type == 'COH2D4':\n"
    "        arrElemData[k*6+0] = 78 #type is cohesive zone\n"
    "        arrElemData[k*6+2] = elems[k].connectivity[0]\n"
    "        arrElemData[k*6+3] = elems[k].connectivity[1]\n"
    "        arrElemData[k*6+4] = elems[k].connectivity[2]\n"
    "        arrElemData[k*6+5] = elems[k].connectivity[3]\n"
    "        nCZS += 1\n"
    "print(\"nElems {0}; nCZS {1}; nTriangles {2}\".format(nElems,nCZS,nTriangles))\n"
    "buf = struct.pack('%si' % len(arrElemData), *arrElemData)\n"
    "f.write(buf)\n";


    s << "arrFrameNodeData = array.array('d', (0 for i in range(0,nNodeValues*4))) #U, V\n"
    "arrFrameNodeLabels = array.array('i', (0 for i in range(0,nNodeValues*2)))\n"
    "arrFrameStressData = array.array('d', (0 for i in range(0,nStressValues*3))) # mises, max principal, min principal\n"
    "arrFrameStressLabelsAndType = array.array('i', (0 for i in range(0,nStressValues*2))) # label; type\n"
    "arrFrameSDEGData = array.array('d', (0 for i in range(0,nSDEGValues))) # only for saving SDEG for CZs\n"
    "arrFrameSDEGLabels = array.array('i', (0 for i in range(0,nSDEGValues))) # element labels for CZs\n"
    "arrFrameStatusData = array.array('d', (0 for i in range(0,nTrianglesInFrame))) # mises, max principal, min principal\n"
    "arrFrameTrisLabels = array.array('i', (0 for i in range(0,nTrianglesInFrame*2))) # label; type\n";


    s << "print(\"nFrames {0}; nNodeValues {1}\".format(nFrames, nNodeValues))\n"

    "for k in range(nFrames):\n"
    "    if(k%2 == 0): print(\"exporting frame {0}\".format(k))\n"
    "    fr = fodb.steps['Step-1'].frames[k]\n"
    "    frTime = fr.frameValue\n"
    "    frUval = fr.fieldOutputs['U'].values\n"
    "    frVval = fr.fieldOutputs['V'].values\n"
    "    frSval = fr.fieldOutputs['S'].values\n"
    "    frSDEGval = fr.fieldOutputs['SDEG'].values\n"
    "    frSTATUSval = fr.fieldOutputs['STATUS'].values\n"
    "    f.write(struct.pack('d', frTime))\n"

    "    # node data per frame\n"
    "    for j in range(nNodeValues):\n"
    "        arrFrameNodeLabels[j*2+0] = frUval[j].nodeLabel\n"
    "        arrFrameNodeData[j*4+0] = frUval[j].dataDouble[0]\n"
    "        arrFrameNodeData[j*4+1] = frUval[j].dataDouble[1]\n"
    "        arrFrameNodeData[j*4+2] = frVval[j].dataDouble[0]\n"
    "        arrFrameNodeData[j*4+3] = frVval[j].dataDouble[1]\n"
    "        if(frUval[j].instance == inst_block):\n"
    "            arrFrameNodeLabels[j*2+1] = 1 # node belongs to the block\n"
    "        else:\n"
    "            arrFrameNodeLabels[j*2+1] = -1\n"

    "    buf1 = struct.pack('%si' % len(arrFrameNodeLabels),*arrFrameNodeLabels)\n"
    "    buf2 = struct.pack('%sd' % len(arrFrameNodeData),*arrFrameNodeData)\n"
    "    f.write(buf1)\n"
    "    f.write(buf2)\n"

    "    # element data per frame\n"
    "    for j in range(nStressValues):\n"
    "        if(frSval[j].baseElementType == 'CPS3'):\n"
    "            arrFrameStressLabelsAndType[j*2+0] = 77\n"
    "            arrFrameStressData[j*3+0] = frSval[j].mises\n"
    "            arrFrameStressData[j*3+1] = frSval[j].maxPrincipal\n"
    "            arrFrameStressData[j*3+2] = frSval[j].minPrincipal\n"
    "        elif(frSval[j].baseElementType == 'COH2D4'):\n"
    "            arrFrameStressLabelsAndType[j*2+0] = 78\n"
    "        else:\n"
    "            arrFrameStressLabelsAndType[j*2+0] = -1\n"
    "        arrFrameStressLabelsAndType[j*2+1] = frSval[j].elementLabel\n"

    "    for j in range(nSDEGValues):\n"
    "        if(frSDEGval[j].baseElementType == 'COH2D4'):\n"
    "            arrFrameSDEGLabels[j] = frSDEGval[j].elementLabel\n"
    "            arrFrameSDEGData[j] = frSDEGval[j].data\n"
    "        else:\n"
    "            arrFrameSDEGLabels[j] = -1\n"

    "    buf1 = struct.pack('%si' % len(arrFrameStressLabelsAndType),*arrFrameStressLabelsAndType)\n"
    "    buf2 = struct.pack('%sd' % len(arrFrameStressData),*arrFrameStressData)\n"
    "    buf3 = struct.pack('%si' % len(arrFrameSDEGLabels),*arrFrameSDEGLabels)\n"
    "    buf4 = struct.pack('%sd' % len(arrFrameSDEGData),*arrFrameSDEGData)\n"
    "    f.write(buf1)\n"
    "    f.write(buf2)\n"
    "    f.write(buf3)\n"
    "    f.write(buf4)\n"

    "    #only triangles - save status vaiable\n"
    "    for j in range(nTrianglesInFrame):\n"
    "        arrFrameTrisLabels[j*2+0] = frSTATUSval[j].elementLabel\n"
    "        if(frSTATUSval[j].baseElementType == 'CPS3'):\n"
    "            arrFrameTrisLabels[j*2+1] = 77\n"
    "        elif(frSTATUSval[j].baseElementType ==  'COH2D4'):\n"
    "            arrFrameTrisLabels[j*2+1] = 78\n"
    "        else:\n"
    "            arrFrameTrisLabels[j*2+1] = -1\n"
    "        arrFrameStatusData[j] = frSTATUSval[j].data\n"
    "    buf1 = struct.pack('%si' % len(arrFrameTrisLabels),*arrFrameTrisLabels)\n"
    "    buf2 = struct.pack('%sd' % len(arrFrameStatusData),*arrFrameStatusData)\n"
    "    f.write(buf1)\n"
    "    f.write(buf2)\n"

    "#close the binary file\n"
    "f.close()\n"
    "print(\"done\")\n";

    s.close();

}
