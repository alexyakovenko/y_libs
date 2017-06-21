
#include "y_neuron.h"


const char * const YNT_FILE_NAME = "nn_3D.ynt"

int main(int argc,char *argv[])
{
register unsigned int _i, _j, _l;
t_ynn *ynn;

//Stage I. Upload NN
//Import ynt
if (!(ynn=import_ynt(YNT_FILE_NAME)))
  { yprintf(YPRINTF_ERROR,"An error occured in import_ynt(), ylib_errno=%s.\n",get_yerrno(ylib_errno)); return EXIT_FAILURE; }
//Compile synapses
if (!(compile_synapses(ynn->size_neurons,ynn->neurons,ynn->size_synapses,ynn->synapses)))
  { yprintf(YPRINTF_ERROR,"An error occured in compile_synapses(), ylib_errno=%s.\n",get_yerrno(ylib_errno)); LABEL_FAILURE_0: free_ynn(ynn), ynn=0x0; return EXIT_FAILURE; }
//Compile neurons
if (!(ynn->size_layers=compile_neurons(ynn->size_neurons,ynn->neurons,ynn->size_synapses,ynn->synapses)))
  { yprintf(YPRINTF_ERROR,"An error occured in compile_neurons(), ylib_errno=%s.\n",get_yerrno(ylib_errno)); goto LABEL_FAILURE_0; }
//Compile layers
if (!(compile_layers(ynn->size_layers,ynn->layers,ynn->size_neurons,ynn->neurons,ynn->size_synapses,ynn->synapses)))
  { yprintf(YPRINTF_ERROR,"An error occured in compile_layers(), ylib_errno=%s.\n",get_yerrno(ylib_errno)); goto LABEL_FAILURE_0; }
else
 {
 ynn->size_output=ynn->layers[ynn->size_layers-1].size_neurons;
 yprintf(YPRINTF_INFO,"ynn is successfully uploaded and compiled.\nBrief statistics is:\n amount of layers is %d\n amount of neurons %d\n amount of synapses %d\nDetails are:\nInput layer has %d neurons\n",ynn->size_layers,ynn->size_neurons,ynn->size_synapses,ynn->size_in);
 for (_i=0; _i<ynn->size_layers; _i++)
   if (ynn->layers[_i].type=='c')
     yprintf(YPRINTF_NOTHING,"Layer %d (CONVOLUTION) consists of %d neurons and %d synapses in %d sublayers\n",_i+1,ynn->layers[_i].size_neurons,ynn->layers[_i].size_synapses,ynn->layars[_i].slayers);
   else 
     yprintf(YPRINTF_NOTHING,"Layer %d (REGULAR) consists of %d neurons and %d synapses\n",_i+1,ynn->layers[_i].size_neurons,ynn->layers[_i].size_synapses);
 }

//Stage II. Upload data
//Import training set
//if (!(import_data(params.train_file_name,&size_train,ynn->size_input,&ids_train,&values_train,&inputs_train,&label_train)))
//  {
//  yprintf(YPRINTF_ERROR,"An error occured during importing train data set, ylib_errno=%s.\n",get_yerrno(ylib_errno));
//  free_ynn(ynn), ynn=0x0;
//  return EXIT_FAILURE;
//  }
//else yprintf(YPRINTF_INFO,"%d samples were imported as training set.\n",pynn->size_train);

//Free (some) memory and exit
free_ynn(ynn); ynn=0x0;

yprintf(YPRINTF_NOTHING,"y_nn is successfully completted!\n");
return EXIT_SUCCESS;
}



