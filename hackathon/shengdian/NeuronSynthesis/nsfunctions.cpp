#include"nsfunctions.h"
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
    vector<int> norder(siz,0);  norder=getNodeOrder(nt);
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
    //debug
//    for(V3DLONG i=0;i<this->listBranch.size();i++)
//    {
//        if(i<10)
//        {
//            cout<<"id:"<<this->listBranch[i].id;
//            cout<<";parent_id:"<<this->listBranch[i].parent_id;
//            cout<<";level:"<<this->listBranch[i].level;
//            cout<<";type:"<<this->listBranch[i].type;
//            cout<<";len:"<<this->listBranch[i].length;
//            cout<<";pathlen:"<<this->listBranch[i].pathLength<<endl;
//        }
//    }
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
NeuronTree branchTree_to_neurontree(const BranchTree& bt)
{
    cout<<"reback the connection"<<endl;
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
NeuronTree tip_branch_pruning(NeuronTree nt, int in_thre)
{
    /*1. get tip, branch and soma nodes;
     * 2. from tip nodes, get tip-branch
     * 3. pruning
     * 4. save
    */
    NeuronTree nt_out;nt_out.listNeuron.clear();nt_out.hashNeuron.clear();
    if(!nt.listNeuron.size())
        return nt_out;
    V3DLONG siz=nt.listNeuron.size();
    /*1. get the index of nt:
                                        * swc_n -> index */
    QHash <int, int>  hashNeuron;hashNeuron.clear();
    V3DLONG somaid=1;
    for (V3DLONG i=0;i<siz;i++)
    {
        hashNeuron.insert(nt.listNeuron[i].n,i);
        if(nt.listNeuron[i].type==1&&nt.listNeuron[i].pn<0)
            somaid=i;
    }
    // 2. get node type: index -> node_type
    /*soma: ntype>=3, branch: ntype=2; tip: ntype=0; internodes: ntype=1*/
    std::vector<int> ntype(siz,0);
    std::vector<int> nkept(siz,1);
    ntype[somaid]=2;
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(s.pn&&hashNeuron.contains(s.pn))
        {
            V3DLONG spn_id=hashNeuron.value(s.pn);
            ntype[spn_id]+=1;
        }
    }
    //3. from tip-nodes
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        std::vector<V3DLONG> nscanned;nscanned.clear();
        if(ntype[i]==0)
        {
            nscanned.push_back(i);
            //distance computing
            double tipbranch_dist=0;
            int ptype=ntype[hashNeuron.value(s.pn)];
            while(ptype)
            {
                V3DLONG spn_id=hashNeuron.value(s.pn);
                NeuronSWC spn=nt.listNeuron[spn_id];
                tipbranch_dist+=sqrt((spn.x-s.x)*(spn.x-s.x)+
                                                        (spn.y-s.y)*(spn.y-s.y)+
                                                        (spn.z-s.z)*(spn.z-s.z));
                if(ptype>1)
                    break;
                nscanned.push_back(spn_id);
                s=spn;ptype=ntype[hashNeuron.value(s.pn)];
            }
            if(tipbranch_dist<=in_thre)
            {
                //remove this tip-branch
                for(V3DLONG n=0;n<nscanned.size();n++)
                    nkept[nscanned[n]]=0;
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
        //debug
//        if(true)
//        {
//            NeuronSWC s = nt.listNeuron[i];
//            if(nkept[i]==0){s.type=1;s.r=10;}
//            nt_out.listNeuron.append(s);
//            nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
//        }
    }

    cout<<"pruning finished"<<endl;
    return nt_out;
}
NeuronTree to_topology_tree(NeuronTree nt)
{
    /*1. get tip, branch and soma nodes;
     * 2. remove internodes
     * 3. reindex
     * 3. save
    */
    NeuronTree nt_out;nt_out.listNeuron.clear();nt_out.hashNeuron.clear();
    if(!nt.listNeuron.size())
        return nt_out;
    V3DLONG siz=nt.listNeuron.size();

    /*1. get the index of nt:
                                        * swc_n -> index */
    QHash <int, int>  hashNeuron;hashNeuron.clear();
    V3DLONG somaid=1;
    for (V3DLONG i=0;i<siz;i++)
    {
        hashNeuron.insert(nt.listNeuron[i].n,i);
        if(nt.listNeuron[i].type==1&&nt.listNeuron[i].pn<0)
            somaid=i;
    }
    // 2. get node type: index -> node_type
    /*soma: ntype>=3, branch: ntype=2; tip: ntype=0; internodes: ntype=1*/
    std::vector<int> ntype(siz,0);
    ntype[somaid]=2;
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(s.pn&&hashNeuron.contains(s.pn))
        {
            V3DLONG spn_id=hashNeuron.value(s.pn);
            ntype[spn_id]+=1;
        }
    }
    //3. remove internode
    //add soma-node first
    NeuronSWC soma_node = nt.listNeuron[somaid];
    nt_out.listNeuron.append(soma_node);
    nt_out.hashNeuron.insert(soma_node.n,nt_out.listNeuron.size()-1);
    //add tip and branch nodes
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        long outpid=s.pn;
        int ptype=ntype[hashNeuron.value(outpid)];
        if(ntype[i]==0||ntype[i]==2)
        {
            //1. for tips to branch-node or soma-node
            //2. for branch-nodes to branch-node or soma-node
            while(ptype==1)
            {
                V3DLONG spn_id=hashNeuron.value(outpid);
                outpid=nt.listNeuron[spn_id].pn;
                ptype=ntype[hashNeuron.value(outpid)];
            }
        }
        else //internode
            continue;
        //update parent id
        s.pn=outpid;
        nt_out.listNeuron.append(s);
        nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
    }
    //pruning
//    NeuronTree nt_out_pruning=tip_branch_pruning(nt_out,0);
    //reindex
//     NeuronTree nt_out_reindex=reindexNT(nt_out);
//    return nt_out_reindex;
    return nt_out;
}
NeuronTree reindexNT(NeuronTree nt)
{
    NeuronTree nt_out_reindex,nt_f;

    QHash <V3DLONG, V3DLONG>  hashNeuron;hashNeuron.clear();
    V3DLONG siz=nt.listNeuron.size();
    V3DLONG somaid=-1;
    for (V3DLONG i=0;i<siz;i++)
    {
        hashNeuron.insert(nt.listNeuron[i].n,i);
        if(nt.listNeuron[i].type==1&&nt.listNeuron[i].pn<0&&somaid<0)
            somaid=i;
    }
    cout<<"soma id "<<somaid<<endl;
    if(somaid<0)
        return nt_out_reindex;

    //only kept nodes that can reach to soma node
    std::vector<int> nkept(siz,-1);    nkept[somaid]=1;
    for(V3DLONG i=0;i<siz;i++){
        NeuronSWC s=nt.listNeuron.at(i);
        if(nkept.at(i)<0){
            if(!hashNeuron.contains(s.pn)) {nkept[i]=0;continue;}
            QList<V3DLONG> scanned_id_list; scanned_id_list.clear();
            scanned_id_list.append(i);
            bool isok=false;
            V3DLONG sp_id=hashNeuron.value(s.pn);
            NeuronSWC sp=nt.listNeuron.at(sp_id);
            while(true){
                if(nkept.at(sp_id)>0) {isok=true; break;}
                scanned_id_list.append(sp_id);
                s=sp;
                if(!hashNeuron.contains(s.pn)) {break;}
                sp_id=hashNeuron.value(s.pn);
                sp=nt.listNeuron.at(sp_id);
            }
            for(V3DLONG sj=0;sj<scanned_id_list.size();sj++)
                nkept[scanned_id_list.at(sj)]=(isok)?1:0;
        }
    }
    for(V3DLONG i=0;i<siz;i++){
        if(nkept.at(i)>0){
            nt_f.listNeuron.append(nt.listNeuron.at(i));
            nt_f.hashNeuron.insert(nt.listNeuron.at(i).n,nt_f.listNeuron.size()-1);
        }
    }

    //reindex
    for(V3DLONG i=0;i<nt_f.listNeuron.size();i++)
    {
        NeuronSWC s = nt_f.listNeuron[i];
        s.pn=(s.pn<0)?s.pn:(nt_f.hashNeuron.value(s.pn)+1);
//        s.n=nt.hashNeuron.value(s.n)+1;
        s.n=i+1;
        nt_out_reindex.listNeuron.append(s);
        nt_out_reindex.hashNeuron.insert(s.n,nt_out_reindex.listNeuron.size()-1);
    }
   return nt_out_reindex;
}
NeuronTree redundancy_bifurcation_pruning(NeuronTree nt,bool not_remove_just_label)
{
    /*1, get node type
     * 2, get three-bifurcation nodes
     * 3, for each three-bifurcation nodes
                *get its child nodes
                * remove the smallest arbors
    */
    NeuronTree nt_out;    V3DLONG siz=nt.listNeuron.size();
    if(!siz)
        return nt_out;
    //get node type
    vector<int> ntype(siz,0);    ntype=getNodeType(nt);
     vector<int> nkept(siz,1);
    //get soma index
    V3DLONG somaid=1;
    vector< vector<V3DLONG> > child_index_list(siz,vector<V3DLONG>());
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(s.pn<0&&s.type==1)
            somaid=i;
        else
        {
            V3DLONG p_index=nt.hashNeuron.value(s.pn);
            child_index_list[p_index].push_back(i);
        }
    }
    //get three bifurcation nodes
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(ntype[i]>2&&i!=somaid)
        {
            cout<<"bifruction node id: "<<s.n<<endl;
            cout<<"child size: "<<child_index_list.at(i).size()<<endl;
            //get all subtrees
            vector<double> child_len_list(child_index_list.size(),0);
            QList<int> kept_child_index_list; kept_child_index_list.clear();
            kept_child_index_list.append(0);kept_child_index_list.append(1);
            for(int c=0;c<child_index_list[i].size();c++)
            {
                V3DLONG this_child_index=child_index_list[i][c];
                NeuronSWC this_sc=nt.listNeuron[this_child_index];
                cout<<"child index "<<c<<":"<<this_sc.n<<endl;
                NeuronTree this_child_subtree=getSubtree(nt,this_sc.n);
                child_len_list[c]=get_nt_len(this_child_subtree);
            }
            //get first kept nodes
            for(int c=1;c<child_index_list[i].size();c++){
                if(child_len_list[c]>child_len_list[kept_child_index_list.at(0)])
                    kept_child_index_list[0]=c;
            }
            //get second kept nodes
            for(int c=0;c<child_index_list[i].size();c++){
                if(child_len_list[c]>child_len_list[kept_child_index_list.at(1)]
                        &&child_len_list[c]<child_len_list[kept_child_index_list.at(0)])
                    kept_child_index_list[1]=c;
            }
            //remove
             for(int c=0;c<child_index_list[i].size();c++){
                 if(c==kept_child_index_list[0]||c==kept_child_index_list[1])
                     continue;
                 V3DLONG remove_child_nt_index=child_index_list[i][c];
                 cout<<"removed arbor len="<<child_len_list[c]<<",id="<<remove_child_nt_index<<endl;
                 NeuronTree remove_tree=getSubtree(nt,nt.listNeuron[remove_child_nt_index].n);
                 nkept[remove_child_nt_index]=0;
                 for(int r=0;r<remove_tree.listNeuron.size();r++)
                 {
                     V3DLONG rc_nt_index=remove_tree.listNeuron[r].n;
                     V3DLONG rc_index=nt.hashNeuron.value(rc_nt_index);
                     nkept[rc_index]=0;
                 }
             }
        }
    }
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(!nkept[i]){
            if(not_remove_just_label)
                s.type=5;
            else
                continue;
        }
        nt_out.listNeuron.append(s);
        nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
    }
    if(!not_remove_just_label)
        nt_out=reindexNT(nt_out);
    return nt_out;
}
NeuronTree three_bifurcation_processing(NeuronTree nt)
{
    //s1. detect three bifurcation points
    //s2. move one of the branch to the parent of bifurcation node
    NeuronTree nt_out;
    V3DLONG siz=nt.listNeuron.size();
    if(!siz)
        return nt_out;
    V3DLONG somaid=1;
    bool another_process=true;
    NeuronTree in_nt;in_nt.deepCopy(nt);
    while(another_process){
        another_process=false;
        //update somaid
        for (V3DLONG i=0;i<siz;i++){
            if(in_nt.listNeuron[i].type==1&&in_nt.listNeuron[i].pn<0){
                somaid=i; break;
            }
        }
        std::vector<int> ntype(siz,0);       ntype[somaid]=2;
        for (V3DLONG i=0;i<siz;i++)
        {
            NeuronSWC s = in_nt.listNeuron[i];
            if(s.pn)
                ntype[in_nt.hashNeuron.value(s.pn)]+=1;
        }
        for (V3DLONG i=0;i<siz;i++)
        {
            NeuronSWC s = in_nt.listNeuron[i];
            if(ntype[i]<=2&&s.pn>0)
            {
                V3DLONG spn_id=in_nt.hashNeuron.value(s.pn);
                if(ntype[spn_id]>2&&spn_id!=somaid)
                {
                    if(ntype[spn_id]>3)
                        another_process=true;
                    cout<<"bifruction node id: "<<s.pn<<endl;
                    V3DLONG spn_id_iter=spn_id;
                    NeuronSWC sp = in_nt.listNeuron[spn_id_iter];
                    while(ntype[spn_id_iter]!=1)
                    {
                        spn_id_iter=in_nt.hashNeuron.value(sp.pn);
                        s.pn=sp.pn;
                        if(spn_id_iter==somaid)
                            break;
                        sp = in_nt.listNeuron[spn_id_iter];
                    }
                    ntype[spn_id]-=1;
                }
            }
            nt_out.listNeuron.append(s);
            nt_out.hashNeuron.insert(s.n,nt_out.listNeuron.size()-1);
        }
        if(another_process)
        {
            in_nt.listNeuron.clear();in_nt.hashNeuron.clear();
            in_nt.deepCopy(nt_out);
            nt_out.listNeuron.clear();nt_out.hashNeuron.clear();
            cout<<"Another iteration"<<endl;
        }
    }
    //check the performance
    return nt_out;
}
bool split_neuron_type(QString inswcpath,QString outpath,int saveESWC)
{
    /*make sure type consistence
     * dendrite,axon,apical_dendrite
     * each type should have soma point
    */
    NeuronTree nt = readSWC_file(inswcpath);
    V3DLONG siz=nt.listNeuron.size();
    if(!siz) {return false;}

    QHash <V3DLONG, V3DLONG>  hashNeuron;hashNeuron.clear();
    V3DLONG somaid=-1;
    for (V3DLONG i=0;i<siz;i++)
    {
        hashNeuron.insert(nt.listNeuron[i].n,i);
        if(nt.listNeuron[i].type==2&&nt.listNeuron[i].pn<0)
        {
            cout<<"Error: multiple -1 nodes"<<endl;
            return false;
        }
        if(nt.listNeuron[i].type==1&&nt.listNeuron[i].pn<0/*&&somaid<0*/){
            if(somaid>=0){cout<<"Error: multiple soma nodes"<<endl;return false;}
            somaid=i;
        }
    }
    if(somaid<0){cout<<"Error: can't find any soma node"<<endl;   return false;}

    //generation
    NeuronSWC soma_node = nt.listNeuron.at(somaid);

    NeuronTree nt_dendrite_whole;
    nt_dendrite_whole.listNeuron.append(soma_node);
    nt_dendrite_whole.hashNeuron.insert(soma_node.n,nt_dendrite_whole.listNeuron.size()-1);

    NeuronTree nt_dendrite;
    nt_dendrite.listNeuron.append(soma_node);
    nt_dendrite.hashNeuron.insert(soma_node.n,nt_dendrite.listNeuron.size()-1);

    NeuronTree nt_apicaldendrite;
    nt_apicaldendrite.listNeuron.append(soma_node);
    nt_apicaldendrite.hashNeuron.insert(soma_node.n,nt_apicaldendrite.listNeuron.size()-1);

    NeuronTree nt_axon;
    nt_axon.listNeuron.append(soma_node);
    nt_axon.hashNeuron.insert(soma_node.n,nt_axon.listNeuron.size()-1);

    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        NeuronSWC sp;
        if(hashNeuron.contains(s.pn))
            sp=nt.listNeuron[hashNeuron.value(s.pn)];
        if(s.type==3||s.type==4)
        {
            if(sp.type!=3&&sp.type!=4&&sp.type!=1)
            {
                s.pn=soma_node.n;
            }
            nt_dendrite_whole.listNeuron.append(s);
            nt_dendrite_whole.hashNeuron.insert(s.n,nt_dendrite_whole.listNeuron.size()-1);
            if(s.type==3)
            {
                nt_dendrite.listNeuron.append(s);
                nt_dendrite.hashNeuron.insert(s.n,nt_dendrite.listNeuron.size()-1);
            }
            else if(s.type==4)
            {
                nt_apicaldendrite.listNeuron.append(s);
                nt_apicaldendrite.hashNeuron.insert(s.n,nt_apicaldendrite.listNeuron.size()-1);
            }
        }
        else if(s.type==2)
        {
            if(sp.type!=2&&sp.type!=1)
            {
//                cout<<"Error in connection"<<endl;
//                return false;
                s.pn=soma_node.n;
//                s.pn=-1; s.type=1;
            }
            nt_axon.listNeuron.append(s);
            nt_axon.hashNeuron.insert(s.n,nt_axon.listNeuron.size()-1);
        }
    }
    //reorder index
    nt_dendrite_whole=reindexNT(nt_dendrite_whole);
    nt_axon=reindexNT(nt_axon);
    nt_dendrite=reindexNT(nt_dendrite);
    nt_apicaldendrite=reindexNT(nt_apicaldendrite);
    //save
    QFileInfo infileinfo(inswcpath);
    QString outfile=(outpath.isEmpty())?(inswcpath):(outpath+"/"+infileinfo.completeBaseName());

    if(saveESWC)
    {
        if(nt_apicaldendrite.listNeuron.size()>1)
        {
            writeESWC_file(outfile+"_wholeDendrite.eswc",nt_dendrite_whole);
            writeESWC_file(outfile+"_apicalDendrite.eswc",nt_apicaldendrite);
        }
        writeESWC_file(outfile+"_dendrite.eswc",nt_dendrite);
        writeESWC_file(outfile+"_axon.eswc",nt_axon);
    }
    else
    {
        if(nt_apicaldendrite.listNeuron.size()>1)
        {
            writeSWC_file(outfile+"_wholeDendrite.swc",nt_dendrite_whole);
            writeSWC_file(outfile+"_apicalDendrite.swc",nt_apicaldendrite);
        }
        writeSWC_file(outfile+"_dendrite.swc",nt_dendrite);
        writeSWC_file(outfile+"_axon.swc",nt_axon);
    }
    return true;
}
NeuronTree smooth_branch_movingAvearage(NeuronTree nt, int smooth_win_size)
{
    /*method: moving avearage
     * 1. get all the branches
     * 2. smooth every nodes in each branch
     * 3.
    */
    NeuronTree nt_smoothed,nt_out;
    nt_out.listNeuron.clear();nt_out.hashNeuron.clear();
    if(!nt.listNeuron.size())
        return nt_out;
    V3DLONG siz=nt.listNeuron.size();
    /*1. get the index of nt:
                                        * swc_n -> index */
    QHash <V3DLONG, V3DLONG>  hashNeuron;hashNeuron.clear();
    V3DLONG somaid=1;
    for (V3DLONG i=0;i<siz;i++)
    {
        hashNeuron.insert(nt.listNeuron[i].n,i);
        if(nt.listNeuron[i].type==1&&nt.listNeuron[i].pn<0)
            somaid=i;
    }
    // 2. get node type: index -> node_type
    /*soma: ntype>=3, branch: ntype=2; tip: ntype=0; internodes: ntype=1*/
    std::vector<int> ntype(siz,0);
    ntype[somaid]=2;
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = nt.listNeuron[i];
        if(s.pn&&hashNeuron.contains(s.pn))
        {
            V3DLONG spn_id=hashNeuron.value(s.pn);
            ntype[spn_id]+=1;
        }
    }
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
