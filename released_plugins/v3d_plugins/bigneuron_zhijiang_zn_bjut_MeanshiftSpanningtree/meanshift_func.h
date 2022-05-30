/* example_func.h
 * This is an example plugin, the header file of plugin functions.
 * 2012-02-10 : by Yinan Wan
 */
 
#ifndef __EXAMPLE_FUNC_H__
#define __EXAMPLE_FUNC_H__

#include <v3d_interface.h>
#include<Gradient.h>
#include "node.h"
struct input_PARA
{
    QString inimg_file;
    V3DLONG channel;
    double prim_distance;//the distance decide to delete covered nodes
    double threshold;//the threshold which used in determining noisy node
    double percentage;//same effect with threshold
};

int meanshift_plugin_vn4(V3DPluginCallback2 &callback, QWidget *parent);
void bf_vn3(QMap<int,Node* > roots,double **weight_result,unsigned char * &img1d,double average_dis,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
void meanshift_vn5(unsigned char * &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r,int iteration);
Node getnode(Node *node);
void printf_path_between_roots(QMultiMap<V3DLONG,QMultiMap<V3DLONG,QList<Node> > > root_path,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
double Map_coordinate(Node current_center,Node target_node,V3DLONG relate_x,V3DLONG relate_y,V3DLONG relate_z);
QMap<V3DLONG,QList<Node> >  build_adj_matrix(QMultiMap<V3DLONG,QMultiMap<V3DLONG,QList<Node> > > path_between_each_roots,QMap<V3DLONG,QList<Node> >  node_con_node);
Node choose_region_vn4(unsigned char * &img1d,struct sphere_model_two_directions sphere_m,Node source_node,Node target_node,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QMap<V3DLONG,QMap<V3DLONG,QList<Node> > > sphere_search(unsigned char * &img1d,QMap<int,Node* > cluster_root,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,double r,int iteration);
QMultiMap<V3DLONG,QMultiMap<V3DLONG,QList<Node> > > spanning_combined_bf_vn3(unsigned char * &img1d,QMap<int,Node* > cluster_root,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,double r,int iteration);
QMultiMap<V3DLONG,QMultiMap<V3DLONG,QList<Node> > > spanning_combined_bf_vn2(unsigned char * &img1d,QMap<int,Node* > cluster_root,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,double r,int iteration);
Node choose_region_vn3(unsigned char * &img1d,struct sphere_model_four_directions sphere_m,Node source_node,Node target_node,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
Node choose_region_vn2(unsigned char * &img1d,struct sphere_model sphere_m,Node source_node,Node target_node,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<Node> change_type_to_QList(QMap<V3DLONG,QMap<V3DLONG,Node> >  roots);
QList<NeuronSWC> fulfill_path_between_roots(QList<NeuronSWC> con_tree,QMap<V3DLONG,QMap<V3DLONG,QList<Node> > > path_between_each_roots,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<NeuronSWC> spanning_tree_algorithm(double** markEdge,QList<Node> seeds);
QList <NeuronSWC > calculate_distance_between_roots(QMap<V3DLONG,QMap<V3DLONG,QList<Node> > > path_between_each_roots,QMap<int,Node* > roots,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
Node choose_region(struct sphere_model sphere_m,Node source_node,Node target_node);
QMultiMap<V3DLONG,QMultiMap<V3DLONG,QList<Node> > > spanning_combined_bf(unsigned char * &img1d,QMap<int,Node* > cluster_root,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,double r,int iteration);
QMap<V3DLONG,Node> merge_cluster_node(QMap<V3DLONG,Node> rootnodes,double distance);

QMap<V3DLONG,QList<Node> >  check_node_conection(QMultiMap<V3DLONG,QMultiMap<V3DLONG,QList<Node> > > path_between_each_roots,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QMap<int,Node*> delete_cluster_node(unsigned char * &img1d,QMap<V3DLONG,QList<Node> >  final_cluster,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG *in_sz,double prim_distance);
 QMap<V3DLONG,QList<Node> >  cluster2newroot(QMap<V3DLONG,QList<V3DLONG> >  covered,QMultiMap<V3DLONG,Node> cluster);
QList<NeuronSWC> spanning_without_bf(QMap<V3DLONG,QMap<V3DLONG,Node> >  roots);

//QMap<V3DLONG,QMap<V3DLONG,Node> >  delete_cluster_node(unsigned char * &img1d,QMap<V3DLONG,QList<Node> >  final_cluster,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG *in_sz,double prim_distance);
 QMap<V3DLONG,QList<Node> >  cluster2newroot(QMap<V3DLONG,QList<V3DLONG> >  covered,QMultiMap<V3DLONG,Node> cluster);
QList<NeuronSWC> spanning_without_bf(QMap<V3DLONG,QMap<V3DLONG,Node> >  roots);

//QMap<V3DLONG,QMap<V3DLONG,Node> > delete_cluster_node(unsigned char * &img1d,QMap<V3DLONG,QList<Node> > final_cluster,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG *in_sz,double prim_distance);
 QMap<V3DLONG,QList<Node> > cluster2newroot(QMap<V3DLONG,QList<V3DLONG> > covered,QMultiMap<V3DLONG,Node> cluster);
QList<NeuronSWC> spanning_without_bf(QMap<V3DLONG,QMap<V3DLONG,Node> > roots);

bool determine_noisy(unsigned char * &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,double threshold,double percentage,int iter_step);
Node choose_next_node(unsigned char * &img1d,Node* pre,V3DLONG next_x,V3DLONG next_y,V3DLONG next_z,QList<Node> choose,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QMap<int,Node*> merge_rootnode_vn2(QMap<int,Node*> &rootnodes,unsigned char * &img1d,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,double distance);
long meanshift_vn6(unsigned char * &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r,int iteration,V3DLONG *in_sz);
long meanshift_vn7(unsigned char * &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r,int iteration,V3DLONG *in_sz);
bool found_final_vn3(unsigned char * &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r);
QList<NeuronSWC> connect_shortest_path_vn5(unsigned char * &img1d,QList<Node> path,Node* begin,Node* end,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG *in_sz);
QList<NeuronSWC> connect_shortest_path_vn4(unsigned char * &img1d,QList<Node> path,Node* begin,Node* end,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<NeuronSWC> smooth_SWC_radius(QList<NeuronSWC> target,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
double fix_radius(QHash<V3DLONG,NeuronSWC> neuron,V3DLONG Parent,QList<V3DLONG> Childs );
QHash<V3DLONG,V3DLONG> Child_Parent_Wan(QList<NeuronSWC> target,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QMultiHash<V3DLONG,V3DLONG> Parent_Child_Wan(QList<NeuronSWC> target,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
void change_type(std::vector<NeuronSWC *> &target,QList<NeuronSWC> source);
QList<NeuronSWC> connect_shortest_path_vn3(unsigned char * &img1d,QList<Node> path,Node* begin,Node* end,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<NeuronSWC> connect_shortest_path_vn2(unsigned char * &img1d,QList<Node> path,Node* begin,Node* end,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<NeuronSWC> connect_shortest_path(unsigned char * &img1d,QList<Node> path,NeuronSWC begin,NeuronSWC end,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
double distance_calculate_vn2(unsigned char * &img1d,QList<Node> path,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<Node> found_path_vn3( QMap<V3DLONG,Node> path_map, Node* temp,Node* temp1,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<Node*> found_path_vn2( QMap<V3DLONG,Node*> path_map, Node* temp,Node* temp1,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
void bf_vn2(QMap<int,Node* > roots,double **weight_result,unsigned char * &img1d,double average_dis,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
double distance_calculate(unsigned char * &img1d,QList<Node*> path,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
void bf(QMap<int,Node* > roots,double **weight_result,unsigned char * &img1d,double average_dis,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
bool contain(QList<Node> *queue,V3DLONG x,V3DLONG y,V3DLONG z);
int meanshift_plugin(V3DPluginCallback2 &callback, QWidget *parent);
int meanshift_plugin_vn2(V3DPluginCallback2 &callback, QWidget *parent);
QList<Node*> meanshift_vn2(unsigned char * &img1d,int *flag,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r,int iteration);
QList<Node*> meanshift_vn3(unsigned char * &img1d,int *flag,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r,int iteration);
void enlarge_radiusof_root_xy(unsigned char * &img1d,QList<Node*> &class_List,Node * &root,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
bool exist(V3DLONG x,V3DLONG y,V3DLONG z,QList<Node*> List);
void meanshift_vn4(unsigned char * &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r,int iteration);
void construct_tree(QMap<int,QList<Node*> > finalclass_node,unsigned char * &img1d,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
int meanshift_plugin_vn4(V3DPluginCallback2 &callback, QWidget *parent,unsigned char* img1d,V3DLONG *in_sz, QString &image_name,bool bmenu,input_PARA &PARA);
int meanshift_plugin(const V3DPluginArgList & input, V3DPluginArgList & output);
void printHelp();
QList <NeuronSWC> DFS_construct(double** markEdge,QList<Node> seeds,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
void merge_spanning_tree(QList<QList <NeuronSWC> > &tree_part);
void bubble_sort(QMap<int,double> &root_r,int n);
void spanning_tree1(QList<Node*> &seeds);
void prepare_write(QList<NeuronSWC> marker_MST_sorted);
QList<Node*> enlarge_radiusof_allnodes_xy(unsigned char * &img1d,QList<Node*> &class_List,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
QList<Node*> trim_nodes(QList<Node*> seed,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
NeuronTree  post_process(NeuronTree nt);
void mark_nodeList(int mark,QList<Node*> &List,int s);
int found_finalnode_mark(V3DLONG x,V3DLONG y,V3DLONG z,QList<Node*> final_List);
bool found_final(unsigned char * &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,V3DLONG r);
int ClassMark(QList<Node *> final_List);
void merge_rootnode(QMap<int,Node*> &rootnodes,unsigned char * &img1d,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
void printSWCByMap_List(QMap<int,QList<Node*> >  List,char * filename);
QMap<V3DLONG,Node*> searchAndConnectByWs(unsigned char* &img1d,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,unsigned int ws);
V3DLONG meanshift(unsigned char* &img1d,V3DLONG x,V3DLONG y,V3DLONG z,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,unsigned int ws);
void printSwcByList(QList<Node*> nodeList,char* path);
void printSwcByMap(QMap<int,Node*> nodeMap, char *path);
void createTree(unsigned char* &img1d,Node* curNode,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z,unsigned int ws);
bool checkConnect(unsigned char* &img1d,Node* begin,Node* end,unsigned int ws,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
void prundUnConnectNode(QMap<V3DLONG,Node*> nodeMap);
QList<Gradient*> gradient(QList<Node*> nodes,unsigned char * &img1d,V3DLONG sz_x,V3DLONG sz_y,V3DLONG sz_z);
#endif
