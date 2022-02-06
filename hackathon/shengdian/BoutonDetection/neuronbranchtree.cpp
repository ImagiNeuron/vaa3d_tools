#include "neuronbranchtree.h"
void BranchUnit::get_features(){
    //length, path length
    //angle
    this->length=dis(this->listNode.at(0),this->listNode.at(this->listNode.size()-1));
    if(this->listNode.size()==2)
        this->pathLength=this->length;
    else
    {
        for(V3DLONG i=0;i<this->listNode.size()-2;i++)
        {
            NeuronSWC snode=this->listNode.at(i);
            NeuronSWC enode=this->listNode.at(i+1);
            this->pathLength+=dis(snode,enode);
        }
    }
}
void BranchTree::get_globalFeatures(){
    if(!this->initialized||this->listBranch.size()==0) {cout<<"Branchtree isn't initialized."<<endl;return;}
    V3DLONG siz=this->listBranch.size();
    vector<int> btype(siz,0);
    btype=this->getBranchType();

    this->total_length=this->total_path_length=0.0;
    this->tip_branches=this->soma_branches=this->max_branch_level=0;
    this->total_branches=siz;

    for (V3DLONG i=0;i<siz;i++)
    {
        BranchUnit bu = this->listBranch.at(i);
        this->total_length+=bu.length;
        this->total_path_length+=bu.pathLength;
        if(btype[i]==0)
            this->tip_branches+=1;
        if(bu.parent_id<0)
            this->soma_branches+=1;
        this->max_branch_level=(bu.level>this->max_branch_level)?bu.level:this->max_branch_level;
    }
}
bool BranchTree::get_enhacedFeatures()
{
     /*1. get branch type
      * 2. get child index of branch
      * 3. get enhanced features
      *      lclength,lcpathLength,rclength,rcpathLength;
             lslength,lspathLength,rslength,rspathLength;
             lstips,rstips;
    */
    if(!this->initialized||this->listBranch.size()==0) {cout<<"Branchtree isn't initialized."<<endl;return false;}
    V3DLONG siz=this->listBranch.size();
    vector<int> btype(siz,0);
    btype=this->getBranchType();

    vector< vector<V3DLONG> > child_index_list(siz,vector<V3DLONG>());
    for (V3DLONG i=0;i<siz;i++)
    {
        BranchUnit bu = this->listBranch.at(i);
        if(bu.parent_id>0)
        {
            V3DLONG p_index=this->hashBranch.value(bu.parent_id);
            child_index_list[p_index].push_back(i);
        }
    }
    //3
    for (V3DLONG i=0;i<siz;i++)
    {
        QList<V3DLONG> subtreeBrlist;
        if(btype[i]>0){
            if(child_index_list.at(i).size()!=2) { cout<<this->listBranch.at(i).id<<" child branch size: "<<child_index_list.at(i).size()<<endl; return false;}
            //left part
            V3DLONG lc_index=child_index_list.at(i).at(0);
            this->listBranch[i].lclength=this->listBranch[lc_index].length;
            this->listBranch[i].lcpathLength=this->listBranch[lc_index].pathLength;
            subtreeBrlist.clear();
            subtreeBrlist=getSubtreeBranches(lc_index);
            for(int sb=0;sb<subtreeBrlist.size();sb++)
            {
                this->listBranch[i].lslength+=this->listBranch[subtreeBrlist.at(sb)].length;
                this->listBranch[i].lspathLength+=this->listBranch[subtreeBrlist.at(sb)].pathLength;
                if(btype[subtreeBrlist.at(sb)]==0)
                    this->listBranch[i].lstips+=1;
            }
            //right part
            V3DLONG rc_index=child_index_list.at(i).at(1);
            this->listBranch[i].rclength=this->listBranch[rc_index].length;
            this->listBranch[i].rcpathLength=this->listBranch[rc_index].pathLength;
            subtreeBrlist.clear();
            subtreeBrlist=getSubtreeBranches(rc_index);
            for(int sb=0;sb<subtreeBrlist.size();sb++)
            {
                this->listBranch[i].rslength+=this->listBranch[subtreeBrlist.at(sb)].length;
                this->listBranch[i].rspathLength+=this->listBranch[subtreeBrlist.at(sb)].pathLength;
                if(btype[subtreeBrlist.at(sb)]==0)
                    this->listBranch[i].rstips+=1;
            }
            //swap left and right according to tip_num
            if(this->listBranch.at(i).lstips<this->listBranch.at(i).rstips)
            {
                std::swap(this->listBranch[i].lclength,this->listBranch[i].rclength);
                std::swap(this->listBranch[i].lcpathLength,this->listBranch[i].rcpathLength);
                std::swap(this->listBranch[i].lslength,this->listBranch[i].rslength);
                std::swap(this->listBranch[i].lspathLength,this->listBranch[i].rspathLength);
                std::swap(this->listBranch[i].lstips,this->listBranch[i].rstips);
            }
        }
    }
    return true;
}
bool BranchTree::init_branch_sequence()
{
    if(!this->listBranch.size()||!this->initialized)
        return false;
    //get tip-branch
    V3DLONG siz=this->listBranch.size();
    vector<int> btype(siz,0); btype=this->getBranchType();
    for(V3DLONG i=0;i<this->listBranch.size();i++)
    {
        BranchUnit bu = this->listBranch[i];
        if(btype[i]==0)
        {
            //start from tip-branch
            BranchSequence brs;
            brs.seqLength+=bu.length;
            brs.seqPathLength+=bu.pathLength;
            brs.listbr.append(i);
            if(bu.parent_id>0)
            {
                V3DLONG bupid=this->hashBranch.value(bu.parent_id);
                BranchUnit bup=this->listBranch[bupid];
                while(true)
                {
                    brs.seqLength+=bup.length;
                    brs.seqPathLength+=bup.pathLength;
                    brs.listbr.append(bupid);
                    if(bup.parent_id<0)
                        break;
                    bupid=this->hashBranch.value(bup.parent_id);
                    bup=this->listBranch[bupid];
                }
            }
            brs.seqType=this->listBranch[i].type;
            brs.seqSize=brs.listbr.size();
            //debug
//            {
//                cout<<"seq size:"<<brs.listbr.size();
//                cout<<";seq length:"<<brs.seqLength;
//                cout<<";seq path length:"<<brs.seqPathLength;
//                cout<<";type:"<<brs.seqType<<endl;
//            }
            this->branchseq.append(brs);
        }
    }
    cout<<"seq size: "<<this->branchseq.size()<<endl;
    return true;
}
bool BranchTree::init(NeuronTree in_nt){
    nt.deepCopy(in_nt);
    //neuron tree to branches
    V3DLONG siz=nt.listNeuron.size();
    if(!siz)
        return false;
    //get node type
    vector<int> ntype(siz,0);    ntype=getNodeType(nt);
    vector<int> norder(siz,0);
    if(!getNodeOrder(nt,norder)) {return false;}
    V3DLONG soma_index=1;
    for (V3DLONG i=0;i<siz;i++){
        if(nt.listNeuron[i].pn<0&&nt.listNeuron[i].type==1)
        {
            soma_index=i;break;
        }
    }
    QList<V3DLONG> br_parent_list;br_parent_list.clear();
    QList<V3DLONG> br_tail_list;br_tail_list.clear();
    for(V3DLONG i=0;i<siz;i++)
    {
        //from tip / branch node to branch / soma node.
        NeuronSWC s = nt.listNeuron[i];
        if(s.pn<0)
            continue;
        if(ntype[i]==0||ntype[i]==2)
        {
            QList<NeuronSWC> bu_nodes; bu_nodes.append(s);
            BranchUnit bru;
            bru.level=norder[i];
            bru.id=this->listBranch.size()+1;
            V3DLONG sp_id=nt.hashNeuron.value(s.pn);
            int ptype=ntype[sp_id];
            if(ptype==2)
            {
                //this branch doesn't have internode.
                NeuronSWC sp=nt.listNeuron[sp_id];
                bu_nodes.append(sp);
            }
            else
            {
                while(true)
                {
                    NeuronSWC sp=nt.listNeuron[sp_id];
                    bu_nodes.append(sp);
                    sp_id=nt.hashNeuron.value(sp.pn);
                    if(soma_index==sp_id)
                    {
                        NeuronSWC sp=nt.listNeuron[sp_id];
                        bu_nodes.append(sp);
                        bru.parent_id=-1;
                        break;
                    }
                    if(ntype[sp_id]==2)
                    {
                        NeuronSWC sp=nt.listNeuron[sp_id];
                        bu_nodes.append(sp);
                        break;
                    }
                    ptype=ntype[sp_id];
                }
            }

            bru.type=bu_nodes.at(0).type;
            //sort and load into Branch struct
            for(V3DLONG b=bu_nodes.size()-1;b>=0;b--)
            {
                NeuronSWC bu_node=bu_nodes[b];
                bu_node.n=bu_nodes.size()-b;
                bu_node.pn=(b==bu_nodes.size()-1)?(-1):(bu_node.n-1);
                bru.listNode.append(bu_node);
                bru.hashNode.insert(bu_node.n,bru.listNode.size()-1);
            }
            //record the parent id of this branch
            br_tail_list.append(i);
            br_parent_list.append(sp_id);
            bru.get_features();
            this->listBranch.append(bru);
            this->hashBranch.insert(bru.id,this->listBranch.size()-1);
        }
    }
    //get branch parent id
    for(V3DLONG i=0;i<this->listBranch.size();i++)
    {
        if(this->listBranch[i].parent_id<0)
            continue;
        //parent-node -> index ->
        for(V3DLONG t=0;t<br_tail_list.size();t++)
            if(br_tail_list.at(t)==br_parent_list.at(i))
                this->listBranch[i].parent_id=this->listBranch.at(t).id;
    }
    cout<<"branch size: "<<this->listBranch.size()<<endl;
    return true;
}
QList<V3DLONG> BranchTree::getSubtreeBranches(V3DLONG inbr_index){
    /*retrieval branch also contains in list*/
    V3DLONG siz=this->listBranch.size();
    QList<V3DLONG> subtreeBrlist; subtreeBrlist.clear();
    if(!siz||!this->initialized)
        return subtreeBrlist;
    vector<int> btype(siz,0); btype=this->getBranchType();
    for(V3DLONG i=0;i<siz;i++)
    {
        //start from tip-branch
        BranchUnit bu = this->listBranch[i];
        if(i==inbr_index) { subtreeBrlist.append(i); continue;}
        if(bu.parent_id<0)
            continue;
        if(btype[i]==0)
        {
            V3DLONG pbr_index=this->hashBranch.value(bu.parent_id);
            BranchUnit pbu=this->listBranch[pbr_index];
            while(true)
            {
                if(pbr_index!=inbr_index)
                {
                    if(pbu.parent_id>0){
                        pbr_index=this->hashBranch.value(pbu.parent_id);
                        pbu=this->listBranch[pbr_index];
                    }
                    else
                        break;
                }
                else {subtreeBrlist.append(i);break;}
            }
        }
    }
    return subtreeBrlist;
}
vector<int> BranchTree::getBranchType()
{
    V3DLONG siz=this->listBranch.size();
    vector<int> btype(siz,0);
    if(!siz||!this->initialized)
        return btype;
    /*soma-branch, interbranch: ntype=2; tip-branch: ntype=0*/
    for (V3DLONG i=0;i<siz;i++)
    {
        BranchUnit bu = this->listBranch[i];
        if(bu.parent_id&&this->hashBranch.contains(bu.parent_id))
        {
            int spn_id=this->hashBranch.value(bu.parent_id);
            btype[spn_id]+=1;
        }
    }
    return btype;
}
bool BranchTree::normalize_branchTree(){
    /*1. get neuron tree length / path_length
     * 2. branch_len /= neuron_len
    */
    V3DLONG siz=this->listBranch.size();
    if(!siz||!this->initialized)
        return false;
    if(this->total_length==0||this->total_path_length==0)
        return false;

    for (V3DLONG i=0;i<siz;i++)
    {
        this->listBranch[i].normalize_len(this->total_length);
        this->listBranch[i].normalize_pathlen(this->total_path_length);
//        if(this->total_branches)
//            this->listBranch[i].normalize_tip(this->total_branches);
    }
    return true;
}
BranchTree readBranchTree_file(const QString& filename)
{
    BranchTree bt;
    QFile brfile(filename);
    if (! brfile.open(QIODevice::ReadOnly | QIODevice::Text))
        return bt;
    else
    {
        int linestate=-1;BranchUnit bu;
        while (! brfile.atEnd())
        {
            char _buf[1000], *buf;
            brfile.readLine(_buf, sizeof(_buf));
            for (buf=_buf; (*buf && *buf==' '); buf++); //skip space
            if (buf[0]=='#')
            {
                 if (buf[1]=='B'&&buf[2]=='R'&&buf[3]=='S'&&buf[4]=='T'&&buf[5]=='A'
                         &&buf[6]=='R'&&buf[7]=='T')
                 {
                     //start a new branch ,state=0
                     linestate+=1;
                 }
                 else if(buf[1]=='#'&&buf[2]=='F'&&buf[3]=='e'&&buf[4]=='a'&&buf[5]=='t'
                         &&buf[6]=='u'&&buf[7]=='r'&&buf[8]=='e'&&buf[9]=='s')
                 {
                     //start to get branch features,state=1
                     linestate+=1;
                 }
                 else if(buf[1]=='#'&&buf[2]=='N'&&buf[3]=='o'&&buf[4]=='d'&&buf[5]=='e'&&buf[6]=='s')
                 {
                     //start to get branch nodes,state=2
                     linestate+=1;
                 }
                 else if(buf[1]=='B'&&buf[2]=='R'&&buf[3]=='E'&&buf[4]=='N'&&buf[5]=='D')
                 {
                     //end this branch
                     linestate=-1;
                     bt.listBranch.append(bu);
                     continue;
                 }
                 else
                     continue;
            }
            if(linestate==0)
            {
                bu.listNode.clear();bu.hashNode.clear();
            }
            else if(linestate==1)
            {
                QStringList qsl = QString(buf).trimmed().split(",");
                if (qsl.size()<6)   continue;
                bu.id=qsl[0].toLong();
                bu.parent_id=qsl[1].toLong();
                bu.type=qsl[2].toInt();
                bu.level=qsl[3].toInt();
                bu.length=qsl[4].toDouble();
                bu.pathLength=qsl[5].toDouble();
            }
            else if(linestate==2)
            {
                NeuronSWC S;
                QStringList qsl = QString(buf).trimmed().split(",");
                if (qsl.size()<7)   continue;
                S.n = qsl[0].toInt();
                S.type = qsl[1].toInt();
                S.x = qsl[2].toFloat();
                S.y = qsl[3].toFloat();
                S.z = qsl[4].toFloat();
                S.r = qsl[5].toFloat();
                S.pn = qsl[6].toInt();
                bu.listNode.append(S);
                bu.hashNode.insert(S.n,bu.listNode.size()-1);
//                bt.nt.listNeuron.append(S);
//                bt.nt.hashNeuron.insert(S.n,bt.nt.listNeuron.size()-1);
            }
            else
                continue;
        }
    }
    bt.initialized=true;
    bt.init_branch_sequence();
    return bt;
}
bool writeBranchTree_file(const QString& filename, const BranchTree& bt,bool enhanced)
{
    /*File Format:
                 * #BRSTART
                 * ##Features or enhanced features
                 * ###id,parent_id,type,level,length,pathLength
                 * ##Nodes
                 * ###n,type,x,y,z,radius,parent
                 * #BREND
    */
    if (filename.isEmpty()||bt.listBranch.size()==0)
        return false;
    QFile tofile(filename);
    if(tofile.exists())
        cout<<"File overwrite to "<<filename.toStdString()<<endl;
    QString confTitle="#This file is used for recording all the branches in a neuron tree (by shengdian).\n";
    QString brstart="#BRSTART\n"; QString brend="#BREND\n";
    QString brfeatures="##Features\n";
    QString brfHead="#Fhead: id,parent_id,type,level,length,pathLength\n";
    if(enhanced)
        brfHead="#Fhead: id,parent_id,type,level,length,pathLength"
                ",lclength,lcpathLength,lslength,lspathLength,lstips"
                ",rclength,rcpathLength,rslength,rspathLength,rstips\n";
    QString brnodes="##Nodes\n";
    QString brnHead="#Nhead:n,type,x,y,z,radius,parent\n";
    if(tofile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        //title
        tofile.write(confTitle.toAscii());
        tofile.write(brfHead.toAscii());
        tofile.write(brnHead.toAscii());
        //inside for each branch
        for(V3DLONG i=0;i<bt.listBranch.size();i++)
        {
            BranchUnit bu = bt.listBranch[i];
            tofile.write(brstart.toAscii());
            tofile.write(brfeatures.toAscii());
            QString brf=QString::number(bu.id);
            brf+=(","+QString::number(bu.parent_id));
            brf+=(","+QString::number(bu.type));
            brf+=(","+QString::number(bu.level));
            brf+=(","+QString::number(bu.length));
            if(enhanced){
                brf+=(","+QString::number(bu.pathLength));
                brf+=(","+QString::number(bu.lclength));
                brf+=(","+QString::number(bu.lcpathLength));
                brf+=(","+QString::number(bu.lslength));
                brf+=(","+QString::number(bu.lspathLength));
                brf+=(","+QString::number(bu.lstips));
                brf+=(","+QString::number(bu.rclength));
                brf+=(","+QString::number(bu.rcpathLength));
                brf+=(","+QString::number(bu.rslength));
                brf+=(","+QString::number(bu.rspathLength));
                brf+=(","+QString::number(bu.rstips)+"\n");
            }
            else
                brf+=(","+QString::number(bu.pathLength)+"\n");
            tofile.write(brf.toAscii());

            tofile.write(brnodes.toAscii());
            for(int j=0;j<bu.listNode.size();j++)
            {
                NeuronSWC s=bu.listNode[j];
                QString brn=(QString::number(s.n));
                brn+=(","+QString::number(s.type));
                brn+=(","+QString::number(s.x));
                brn+=(","+QString::number(s.y));
                brn+=(","+QString::number(s.z));
                brn+=(","+QString::number(s.r));
                brn+=(","+QString::number(s.parent)+"\n");
                tofile.write(brn.toAscii());
            }
            tofile.write(brend.toAscii());
        }
        tofile.close();
        return true;
    }
    return false;
}
bool writeBranchSequence_file(const QString& filename, const BranchTree& bt,bool enhanced)
{
    /*File Format:
     *  (from soma to tip branch)
                 * #BRSSTART
                 * ##(soma-branch) id,parent_id,type,level,length,pathLength
                 * ...
                 * ##(tip-branch) id,parent_id,type,level,length,pathLength
                 * #BRSEND
    */
    if (filename.isEmpty()||bt.listBranch.size()==0)
        return false;
    QFile tofile(filename);
    if(tofile.exists())
        cout<<"File overwrite to "<<filename.toStdString()<<endl;
    QString confTitle="#This file is used for recording all the branch sequences in a neuron tree (by shengdian).\n";
    QString brsstart="#BRSSTART\n"; QString brsend="#BRSEND\n";
    QString brfHead="#Fhead: id,parent_id,type,level,length,pathLength\n";
    if(enhanced)
        brfHead="#Fhead: id,parent_id,type,level,length,pathLength"
                ",lclength,lcpathLength,lslength,lspathLength,lstips"
                ",rclength,rcpathLength,rslength,rspathLength,rstips\n";
    if(tofile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        //title
        tofile.write(confTitle.toAscii());
        tofile.write(brfHead.toAscii());
        //inside for each branch
        for(V3DLONG i=0;i<bt.branchseq.size();i++)
        {
            tofile.write(brsstart.toAscii());
            BranchSequence brs=bt.branchseq.at(i);
            for(int j=brs.listbr.size()-1;j>=0;j--)
            {
                V3DLONG bid=brs.listbr.at(j);
                BranchUnit bu = bt.listBranch[bid];
                QString brf=QString::number(bu.id);
                brf+=(","+QString::number(bu.parent_id));
                brf+=(","+QString::number(bu.type));
                brf+=(","+QString::number(bu.level));
                brf+=(","+QString::number(bu.length));
                if(enhanced){
                    brf+=(","+QString::number(bu.pathLength));
                    brf+=(","+QString::number(bu.lclength));
                    brf+=(","+QString::number(bu.lcpathLength));
                    brf+=(","+QString::number(bu.lslength));
                    brf+=(","+QString::number(bu.lspathLength));
                    brf+=(","+QString::number(bu.lstips));
                    brf+=(","+QString::number(bu.rclength));
                    brf+=(","+QString::number(bu.rcpathLength));
                    brf+=(","+QString::number(bu.rslength));
                    brf+=(","+QString::number(bu.rspathLength));
                    brf+=(","+QString::number(bu.rstips)+"\n");
                }
                else
                    brf+=(","+QString::number(bu.pathLength)+"\n");
                tofile.write(brf.toAscii());
            }
            tofile.write(brsend.toAscii());
        }
        tofile.close();
        return true;
    }
    return false;
}
bool getNodeOrder(NeuronTree nt,vector<int> & norder)
{
    /*soma order=0
     * Workflow
     * 1. get node type;
     * 2. from one node to soma,count how many branch nodes will be scanned.
     * 3.out
     * PS: neuron tree must have only one soma node
    */
    V3DLONG siz=nt.listNeuron.size(); if(!siz) { return false;}
    vector<int> ntype(siz,0);    ntype=getNodeType(nt);

    QHash <V3DLONG, V3DLONG>  hashNeuron; hashNeuron.clear();
    V3DLONG somaid=get_soma(nt);
    if(somaid<0){return false;}
    for (V3DLONG i=0;i<siz;i++)
    {
        hashNeuron.insert(nt.listNeuron[i].n,i);
        if(ntype.at(i)>2&&somaid!=i)
            return false;
    }
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(somaid==i)
            continue;
        NeuronSWC s_iter=s;
        long pIndex=hashNeuron.value(s_iter.pn);
        int ptype=ntype[pIndex];
        /*for all the nodes except soma_node*/
        while(ptype<3)
        {
            if(ptype==2)
                norder[i]+=1;
            s_iter=nt.listNeuron[pIndex];
            pIndex=hashNeuron.value(s_iter.pn);
            ptype=ntype[pIndex];
        }
        norder[i]+=1;
    }
    return true;
}
std::vector<int> getNodeType(NeuronTree nt)
{
    /*soma: ntype>=3, branch: ntype=2; tip: ntype=0; internodes: ntype=1
    PS: not have to be a single tree */
    V3DLONG siz=nt.listNeuron.size();
    std::vector<int> ntype(siz,0);    if(!siz) {return ntype;}
//    cout<<"size="<<ntype.size()<<endl;
    V3DLONG somaid=get_soma(nt);
    if(somaid>=0)
        ntype[somaid]=2;
    else
        cout<<"no soma node"<<endl;
    /*1. get the index of nt:     * swc_n -> index */
    QHash <V3DLONG, V3DLONG>  hashNeuron;
    for (V3DLONG i=0;i<siz;i++)
        hashNeuron.insert(nt.listNeuron[i].n,i);
    // 2. get node type: index -> node_type
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(s.pn>=0&&hashNeuron.contains(s.pn))
        {
            V3DLONG spn_id=hashNeuron.value(s.pn);
            ntype[spn_id]+=1;
        }
    }

//    cout<<"size="<<ntype.size()<<endl;
    return ntype;
}
NeuronTree reindexNT(NeuronTree nt)
{
    /*if parent node not exist, will set to -1*/
    NeuronTree nt_out_reindex;
    for(V3DLONG i=0;i<nt.listNeuron.size();i++)
    {
        NeuronSWC s = nt.listNeuron.at(i);
        s.pn=(s.pn<0||!nt.hashNeuron.contains(s.pn))?(-1):(nt.hashNeuron.value(s.pn)+1);
        s.n=i+1;
        nt_out_reindex.listNeuron.append(s);
        nt_out_reindex.hashNeuron.insert(s.n,nt_out_reindex.listNeuron.size()-1);
    }
   return nt_out_reindex;
}
double getNT_len(NeuronTree nt,float *res)
{
    if(res[0]<=0||res[1]<=0||res[2]<=0){
        res[0]=res[1]=res[2]=1.0;
    }
    double out_len=0.0;
    V3DLONG siz=nt.listNeuron.size();
    if(!siz)
        return out_len;
     for (V3DLONG i=0;i<siz;i++)
     {
          NeuronSWC s = nt.listNeuron[i];
          if(s.pn>0&&nt.hashNeuron.contains(s.pn))
          {
              V3DLONG spid=nt.hashNeuron.value(s.pn);
              NeuronSWC sp=nt.listNeuron[spid];
              out_len+=sqrt
                      (res[0]*res[0]*(s.x-sp.x)*(s.x-sp.x)+
                      res[1]*res[1]*(s.y-sp.y)*(s.y-sp.y)+
                      res[2]*res[2]* (s.z-sp.z)*(s.z-sp.z));
          }
     }
    return out_len;
}
/*swc processing*/
V3DLONG get_soma(NeuronTree & nt,bool connect){
    V3DLONG niz=nt.listNeuron.size();
    V3DLONG somaid=-1;
    if(connect){
        for(V3DLONG i=0;i<niz;i++){
            NeuronSWC s=nt.listNeuron.at(i);
            if(s.pn<0&&s.type!=1){
                cout<<"---------------Attempt to process multiple -1 nodes-----------------------"<<endl;
                //find the node with same coordinates
                for(V3DLONG j=0;j<niz;j++){
                    NeuronSWC sj=nt.listNeuron.at(j);
                    if(i!=j&&s.x==sj.x&&s.y==sj.y&&s.z==sj.z)
                    {
                        nt.listNeuron[i].pn=sj.n;
                    }
                }
            }
        }
    }
    for(V3DLONG i=0;i<niz;i++){
        NeuronSWC s=nt.listNeuron.at(i);
        if(s.pn<0){
            if(s.type==1){
                if(somaid>0)
                {
                    cout<<"---------------Error: multiple soma nodes!!!-----------------------"<<endl;
                    return -1;
                }else
                    somaid=i;
            }
            else{
                cout<<"-------------- multiple -1 nodes!!!-----------------------"<<endl;
                return -1;
            }
        }
    }
    return somaid;
}
bool loop_checking(NeuronTree nt){
    V3DLONG siz=nt.listNeuron.size();
    if(!siz) {return false;}
    QHash <V3DLONG, V3DLONG>  hashNeuron;
    for (V3DLONG i=0;i<siz;i++)
        hashNeuron.insert(nt.listNeuron.at(i).n,i);

    QVector<V3DLONG> scanned(siz,0);
    int loop_count=0;
    for (V3DLONG i=0;i<siz;i++){
        NeuronSWC s = nt.listNeuron.at(i);
        if(scanned.at(i)>0)
            continue;
        scanned[i]=1;
        if(s.pn>0&&hashNeuron.contains(s.pn)){
            QList<V3DLONG> snodes;
            snodes.clear();
            snodes.append(i);
            V3DLONG pid=hashNeuron.value(s.pn);
            NeuronSWC sp=nt.listNeuron.at(pid);
            while(true){
                scanned[pid]=1;
                snodes.append(pid);
                if(s.n==sp.n)
                {
                    loop_count++;
                    for(int s=0;s<snodes.size();s++)
                        scanned[snodes.at(s)]=loop_count+1;
                    cout<<"Loop at node "<<s.n<<endl;
                    break;
                }
                if(sp.pn<=0 || !hashNeuron.contains(sp.pn))
                    break;
                pid=hashNeuron.value(sp.pn);
                sp=nt.listNeuron.at(pid);
            }
        }
    }
    if(loop_count>0){
        cout<<"Total loop="<<loop_count<<endl;
        return true;
    }
    return false;
}
bool three_bifurcation_processing(NeuronTree& in_nt)
{
    //s1. detect three bifurcation points
    //s2. move one of the branch to the parent of bifurcation node
    //process bifurcation nodes one by one

    V3DLONG siz=in_nt.listNeuron.size();
    if(!siz) {return false;}

    V3DLONG somaid=get_soma(in_nt,true);
    if(somaid<0){cout<<"Soma error"<<endl; return false;}
    std::vector<int> ntype(siz,0);    ntype=getNodeType(in_nt);

    //s1
    QList<V3DLONG> bifur_idlist; bifur_idlist.clear();
    for (V3DLONG i=0;i<siz;i++)
        if(somaid!=i&&ntype.at(i)>2)
            bifur_idlist.append(i);

    cout<<"#bifurcation="<<bifur_idlist.size()<<endl;
    //s2
    for(int b=0;b<bifur_idlist.size();b++){
        V3DLONG bifur_id=bifur_idlist.at(b);
        bool processed=false;
        for (V3DLONG i=0;i<siz;i++){
            NeuronSWC s = in_nt.listNeuron.at(i);
            if(s.pn>0&&in_nt.hashNeuron.contains(s.pn)){
                 V3DLONG spn_id=in_nt.hashNeuron.value(s.pn);
                 if(spn_id==bifur_id){
                     //find the child node
                     cout<<"Node: "<<s.pn<<",#child="<<ntype.at(spn_id)<<endl;
                     NeuronSWC sp = in_nt.listNeuron[bifur_id];
                     while(true){
                         s=sp;
                         if(s.pn<0||!in_nt.hashNeuron.contains(s.pn)){
                             processed=false;
                             break;
                         }
                         spn_id=in_nt.hashNeuron.value(s.pn);
                         sp = in_nt.listNeuron.at(spn_id);
                         if(ntype.at(spn_id)==1&&spn_id!=somaid){
                             in_nt.listNeuron[i].pn=sp.n;
                             cout<<"set "<<in_nt.listNeuron[i].n<<" parent to "<<sp.n<<endl;
                             processed=true;
                             ntype[spn_id]+=1;
                             break;
                         }
                     }
                 }
            }
            if(processed)
                break;
        }
        if(!processed){cout<<"process fail"<<endl; return false;}
    }
    //check
    ntype=getNodeType(in_nt);
    for (V3DLONG i=0;i<siz;i++)
        if(somaid!=i&&ntype.at(i)>2){
            cout<<ntype.at(i)<<" childs,";
            cout<<"bifurcation id="<<i<<endl;
            return /*false*/three_bifurcation_processing(in_nt);
        }
    return true;
}
NeuronTree tip_branch_pruning(NeuronTree nt, float in_thre)
{
    /*1. get tip, branch and soma nodes;
     * 2. from tip nodes, get tip-branch
     * 3. pruning
     * 4. save
    */
    NeuronTree nt_out;
    V3DLONG siz=nt.listNeuron.size(); if(!siz) {return nt_out;}
    QHash <V3DLONG, V3DLONG>  hashNeuron;hashNeuron.clear();
    for (V3DLONG i=0;i<siz;i++)
        hashNeuron.insert(nt.listNeuron[i].n,i);

    std::vector<int> ntype(siz,0); ntype=getNodeType(nt);
    std::vector<int> nkept(siz,1);

    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron.at(i);
        std::vector<V3DLONG> nscanned;
        nscanned.clear();
        if(ntype.at(i)==0&&s.pn&&hashNeuron.contains(s.pn))
        {
            nscanned.push_back(i);
            //distance computing
            V3DLONG spn_id=hashNeuron.value(s.pn);
            NeuronSWC spn=nt.listNeuron.at(spn_id);
            double tipbranch_dist=0.0;
            while(true)
            {
                tipbranch_dist+=sqrt((spn.x-s.x)*(spn.x-s.x)+
                                                        (spn.y-s.y)*(spn.y-s.y)+
                                                        (spn.z-s.z)*(spn.z-s.z));
                if(ntype.at(spn_id)>1) {break;}
                nscanned.push_back(spn_id);
                s=spn;
                if(s.pn<0||!hashNeuron.contains(s.pn)){
                    nscanned.clear(); tipbranch_dist+=(in_thre+1);break;
                }
                spn_id=hashNeuron.value(s.pn);
                spn=nt.listNeuron.at(spn_id);
            }
            if(tipbranch_dist<=in_thre)
            {
                //remove this tip-branch
                cout<<"tip branch len="<<tipbranch_dist;
                cout<<",small tip-branch: "<<nt.listNeuron.at(i).n<<endl;

                for(V3DLONG n=0;n<nscanned.size();n++)
                    nkept[nscanned.at(n)]=0;
            }
        }
        else //internode
            continue;
    }
    //save
    for (V3DLONG i=0;i<siz;i++)
    {
        if(nkept[i]!=0)
        {
            NeuronSWC s = nt.listNeuron[i];
            nt_out.listNeuron.append(s);
            nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
        }
    }
    cout<<"# pruning nodes="<<nt.listNeuron.size()-nt_out.listNeuron.size()<<endl;
    return nt_out;
}
NeuronTree node_interpolation(NeuronTree nt,int Min_Interpolation_Pixels,bool sort_index){
    cout<<"linear interpolation of neuron tree"<<endl;
    V3DLONG siz = nt.listNeuron.size(); if(!siz) return nt;

    V3DLONG max_index=siz;
    QList<NeuronSWC> listNeuron =  nt.listNeuron;
    QHash <V3DLONG, V3DLONG>  hashNeuron;
    for(V3DLONG i=0;i<siz;i++){
        hashNeuron.insert(listNeuron[i].n,i);
        max_index=MAX(max_index,listNeuron.at(i).n);
    }
    //step2
    QList <NeuronSWC> nt_out_listNeuron;
    V3DLONG new_node_count=0;
    for (V3DLONG i=0;i<listNeuron.size();i++)
    {
        NeuronSWC s = listNeuron.at(i);
        if(s.parent>0&&hashNeuron.contains(s.parent))
        {
            V3DLONG pid=hashNeuron.value(s.parent);
            NeuronSWC sp=listNeuron.at(pid);
            double cp_dist=dis(s,sp);
            double Min_Interpolation_Pixels_dist=sqrt(long(Min_Interpolation_Pixels*Min_Interpolation_Pixels));
            int interpolate_times=int(cp_dist/Min_Interpolation_Pixels_dist);

            if(interpolate_times==1)
            {
                NeuronSWC s_interpolated=s;
                s_interpolated.n=max_index+1;
                //one node at the center of p and sp
                s_interpolated.x=(sp.x+s.x)/float(interpolate_times+1);
                s_interpolated.y=(sp.y+s.y)/float(interpolate_times+1);
                s_interpolated.z=(sp.z+s.z)/float(interpolate_times+1);
                s.parent=s_interpolated.n;
                max_index++;
                nt_out_listNeuron.append(s_interpolated);
            }
            else if(interpolate_times>1)
            {
                //interpolate list of nodes
                float x_Interpolation_dis=(sp.x-s.x)/float(interpolate_times);
                float y_Interpolation_dis=(sp.y-s.y)/float(interpolate_times);
                float z_Interpolation_dis=(sp.z-s.z)/float(interpolate_times);

                NeuronSWC s_interpolated_start=s;
                long spid=s.pn;
                for(int ti=1;ti<=interpolate_times;ti++)
                {
                    NeuronSWC s_interpolated=s_interpolated_start;
                    s_interpolated.n=max_index+1;
                    s_interpolated.x=s_interpolated_start.x+x_Interpolation_dis;
                    s_interpolated.y=s_interpolated_start.y+y_Interpolation_dis;
                    s_interpolated.z=s_interpolated_start.z+z_Interpolation_dis;
                    s_interpolated.parent=max_index+2;
                    if(ti==interpolate_times)
                        s_interpolated.parent=spid;
                    else if(ti==1)
                        s.parent=s_interpolated.n;
                    nt_out_listNeuron.append(s_interpolated);
                    max_index++;
                    s_interpolated_start=s_interpolated;
                }
            }
        }
        nt_out_listNeuron.append(s);
    }
    cout<<"finished the interpolation"<<endl;
    cout<<"from size: "<<siz<<" to "<<nt_out_listNeuron.size()<<endl;
    //step3: re_sort index of nt_out

    NeuronTree nt_out;
    for(V3DLONG i=0;i<nt_out_listNeuron.size();i++)
    {
        NeuronSWC s = nt_out_listNeuron[i];
        nt_out.listNeuron.append(s);
        nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
    }
    if(sort_index)
        return reindexNT(nt_out);
    else
        return nt_out;
}
NeuronTree internode_pruning(NeuronTree nt,float pruning_dist,bool profiled){
    /*1. if pruning_dist<=0, keep branch nodes, soma node, tip nodes
     * 2. if pruning_dist>0, pruning the internode distance below pruning_dist
    */
    V3DLONG siz=nt.listNeuron.size();
    NeuronTree out;    if(siz<=0){return out;}
    QHash <V3DLONG, V3DLONG>  hashNeuron;
    for(V3DLONG i=0;i<siz;i++){
        NeuronSWC s = nt.listNeuron.at(i);
        hashNeuron.insert(s.n,i);
    }
    vector<int> ntype=getNodeType(nt);
    vector<int> nprocessed(siz,0);

    for(V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron.at(i);
        if(nprocessed.at(i)==1)
            continue;
        if(s.parent>0&&hashNeuron.contains(s.parent)
                &&ntype.at(hashNeuron.value(s.parent))==1)
        {
            V3DLONG pid=hashNeuron.value(s.parent);
            NeuronSWC sp=nt.listNeuron.at(pid);
            if(sp.parent>0&&hashNeuron.contains(sp.parent))
            {
                V3DLONG gpid=hashNeuron.value(sp.parent);
                NeuronSWC sgp=nt.listNeuron.at(gpid);
                double ssp_dist=dis(sp,sgp);
                QList<V3DLONG> gplist;gplist.clear();
                while(true)
                {
                    if(ssp_dist>pruning_dist||ntype.at(gpid)>1){
                        break;
                    }
                    gplist.append(gpid);
                    if(sgp.parent<0||!hashNeuron.contains(sgp.parent)
                            ||ntype.at(hashNeuron.value(sgp.parent))>1)
                        break;
                    gpid=hashNeuron.value(sgp.parent);
                    sgp=nt.listNeuron.at(gpid);
                    ssp_dist=dis(sp,sgp);
                }

                if(!gplist.size())
                    continue;
                V3DLONG kept_id=pid; float max_r=sp.r;
                for(int p=0;p<gplist.size();p++)
                {
                    V3DLONG gpid=gplist.at(p);
                    NeuronSWC sgp=nt.listNeuron.at(gpid);
                    nprocessed[gpid]=1;
                    if(profiled&&sgp.r>max_r)
                        kept_id=gpid;
                }
                if(kept_id==pid)
                {
                    //only keep parent node
                    V3DLONG gpid=gplist.at(gplist.size()-1);
                    NeuronSWC sgp=nt.listNeuron.at(gpid);
                    nt.listNeuron[kept_id].parent=sgp.parent;
                }
                else
                {
                    //keep one of the grand parent node
//                    nt.listNeuron[i].parent=sp.parent;
                    nprocessed[pid]=1;
                    //s connect to this node
                    gpid=kept_id;
                    NeuronSWC sgp=nt.listNeuron.at(gpid);
                    nt.listNeuron[i].parent=sgp.n;
                    //this node to end node
                    V3DLONG egpid=gplist.at(gplist.size()-1);
                    NeuronSWC segp=nt.listNeuron.at(egpid);
                    nt.listNeuron[gpid].parent=segp.parent;
                    nprocessed[gpid]=0;
                }
            }
        }
    }
    //for tip nodes
    for(V3DLONG i=0;i<siz;i++){
        NeuronSWC s = nt.listNeuron.at(i);
        if(nprocessed[i]==1)
            continue;
        if(s.parent>0&&hashNeuron.contains(s.parent)){
            V3DLONG pid=hashNeuron.value(s.parent);
            NeuronSWC sp=nt.listNeuron.at(pid);
            if(sp.parent>0&&hashNeuron.contains(sp.parent)
                    &&ntype.at(i)==0){
                double ssp_dist=dis(sp,s);
                if(ssp_dist<pruning_dist)
                {
                    if(profiled&&s.r>=sp.r
                            &&ntype.at(pid)==1
                            &&ntype.at(hashNeuron.value(sp.parent))==1){
                        //consider radius feature: keep nodes with bigger radius
                        //keep tip node
                        nt.listNeuron[i].parent=sp.parent;
                        nprocessed[pid]=1;
                    }
                    else//keep parent node
                        nprocessed[i]=1;
                }
            }
        }
    }

    for(V3DLONG i=0;i<siz;i++){
        NeuronSWC s = nt.listNeuron.at(i);
        if(nprocessed.at(i)==0){
            out.listNeuron.append(s);
            out.hashNeuron.insert(s.n,out.listNeuron.size()-1);
        }
    }
    cout<<"pruning size="<<(nt.listNeuron.size()-out.listNeuron.size())<<endl;
    return reindexNT(out);
}
NeuronTree duplicated_tip_branch_pruning(NeuronTree nt,float dist_thre){
    NeuronTree nt_out;
    V3DLONG siz=nt.listNeuron.size(); if(!siz) {return nt_out;}
    QHash <V3DLONG, V3DLONG>  hashNeuron;hashNeuron.clear();
    for (V3DLONG i=0;i<siz;i++)
        hashNeuron.insert(nt.listNeuron[i].n,i);

    std::vector<int> ntype(siz,0); ntype=getNodeType(nt);
    std::vector<int> nkept(siz,1);

    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron.at(i);
        std::vector<V3DLONG> nscanned;
        nscanned.clear();
        if(ntype.at(i)==0&&s.pn&&hashNeuron.contains(s.pn))
        {
            nscanned.push_back(i);
            //distance computing
            V3DLONG spn_id=hashNeuron.value(s.pn);
            NeuronSWC spn=nt.listNeuron.at(spn_id);
            double tipbranch_dist=0.0;
            while(true)
            {
                tipbranch_dist+=sqrt((spn.x-s.x)*(spn.x-s.x)+
                                                        (spn.y-s.y)*(spn.y-s.y)+
                                                        (spn.z-s.z)*(spn.z-s.z));
                if(ntype.at(spn_id)>1) {break;}
                nscanned.push_back(spn_id);
                s=spn;
                if(s.pn<0||!hashNeuron.contains(s.pn)){
                    nscanned.clear(); tipbranch_dist+=(dist_thre+1);break;
                }
                spn_id=hashNeuron.value(s.pn);
                spn=nt.listNeuron.at(spn_id);
            }
            if(tipbranch_dist<=dist_thre)
            {
                //remove this tip-branch
                cout<<"tip branch len="<<tipbranch_dist;
                cout<<",small tip-branch: "<<nt.listNeuron.at(i).n<<endl;
                bool duplicated=true;
                for(V3DLONG n=0;n<nscanned.size();n++)
                {
                    NeuronSWC cs=nt.listNeuron.at(nscanned.at(n));
                    for(V3DLONG ii=0;ii<siz;ii++){
                        NeuronSWC si=nt.listNeuron.at(ii);
                        if(nscanned.at(n)!=ii&&
                                cs.x==si.x&&
                                cs.y==si.y&&
                                cs.z==si.z){
                            duplicated=true;
                            break;
                        }
                    }
                    if(!duplicated)
                        break;
                }
                if(!duplicated)
                    nscanned.clear();
                else
                    for(V3DLONG n=0;n<nscanned.size();n++)
                        nkept[nscanned.at(n)]=0;
            }
        }
        else //internode
            continue;
    }
    //save
    for (V3DLONG i=0;i<siz;i++)
    {
        if(nkept[i]!=0)
        {
            NeuronSWC s = nt.listNeuron[i];
            nt_out.listNeuron.append(s);
            nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
        }
    }
    cout<<"#dup branch pruning nodes="<<nt.listNeuron.size()-nt_out.listNeuron.size()<<endl;
    return nt_out;
}
NeuronTree smooth_branch_movingAvearage(NeuronTree nt, int smooth_win_size)
{
    /*method: moving avearage
     * 1. get all the branches
     * 2. smooth every nodes in each branch
     * 3.
    */
    NeuronTree nt_smoothed,nt_out;
    V3DLONG siz=nt.listNeuron.size();
    if(!siz)
        return nt_out;
    QHash <V3DLONG, V3DLONG>  hashNeuron;hashNeuron.clear();
    for (V3DLONG i=0;i<siz;i++)
        hashNeuron.insert(nt.listNeuron[i].n,i);

    // 2. get node type: index -> node_type
    /*soma: ntype>=3, branch: ntype=2; tip: ntype=0; internodes: ntype=1*/
    std::vector<int> ntype(siz,0); ntype=getNodeType(nt);
    //get child id: n -> child_n; only for internodes
    QHash <V3DLONG, V3DLONG>  hashChild;hashChild.clear();
    for(V3DLONG i=0;i<siz;i++)
    {
        if(ntype[i]==1)
        {
            V3DLONG this_node_n=nt.listNeuron[i].n;
            for(V3DLONG j=0;j<siz;j++)
            {
                if(i==j)
                    continue;
                if(this_node_n==nt.listNeuron[j].pn)
                {hashChild.insert(i,j); break;}
            }
        }
    }
    //3. start from one key points and end at one key points
    int half_win_size=(smooth_win_size-1)/2;
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        // make no changes to tip, branch and soma nodes
        if(ntype[i]==0||ntype[i]>=2)
        {
            nt_smoothed.listNeuron.append(s);
            nt_smoothed.hashNeuron.insert(s.n,nt_smoothed.listNeuron.size()-1);
            continue;
        }
        int this_smooth_win_size=1;
        //for children of this node
        V3DLONG child_index=hashChild.value(i);
        for(int c=0;c<half_win_size;c++)
        {
            NeuronSWC schild=nt.listNeuron[child_index];
            s.x+=schild.x;
            s.y+=schild.y;
            s.z+=schild.z;
            this_smooth_win_size++;
            if(ntype[child_index]==0||ntype[child_index]>=2)
                break;
            child_index=hashChild.value(child_index);
        }
        //for parents of this node
         V3DLONG parent_index=hashNeuron.value(s.pn);
        for(int p=0;p<half_win_size;p++)
        {
            NeuronSWC sparent=nt.listNeuron[parent_index];
            s.x+=sparent.x;
            s.y+=sparent.y;
            s.z+=sparent.z;
            this_smooth_win_size++;
            if(ntype[parent_index]==0||ntype[parent_index]>=2)
                break;
            parent_index=hashNeuron.value(sparent.pn);
        }
        // avearage
        s.x/=this_smooth_win_size;
        s.y/=this_smooth_win_size;
        s.z/=this_smooth_win_size;
        nt_smoothed.listNeuron.append(s);
        nt_smoothed.hashNeuron.insert(s.n,nt_smoothed.listNeuron.size()-1);
    }
    //remove very near nodes
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt_smoothed.listNeuron[i];
        if(ntype[i]==0||ntype[i]>=2)
            continue;

        NeuronSWC spn=nt_smoothed.listNeuron[hashNeuron.value(s.pn)];
        double to_parent_dist=0;
        to_parent_dist=sqrt((spn.x-s.x)*(spn.x-s.x)+
                                                (spn.y-s.y)*(spn.y-s.y)+
                                                (spn.z-s.z)*(spn.z-s.z));
        if(to_parent_dist<1)
        {
            //remove this node
            //get child id
            V3DLONG child_index=hashChild.value(i);
            nt_smoothed.listNeuron[child_index].pn=s.pn;
            nt_smoothed.listNeuron[i].pn=0;
            continue;
        }
    }
    for (V3DLONG i=0;i<siz;i++)
    {
        if(nt_smoothed.listNeuron[i].pn==0)
            continue;
        NeuronSWC s = nt_smoothed.listNeuron[i];
        nt_out.listNeuron.append(s);
        nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
    }
    cout<<"Remove duplicated node size: "<<(siz-nt_out.listNeuron.size())<<endl;
    qDebug()<<"Finished smoothing process";
    return nt_out;
}
