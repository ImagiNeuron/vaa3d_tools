/*
Fast Artificial Neural Network Library (fann)
Copyright (C) 2003 Steffen Nissen (lukesky@diku.dk)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/  
	
/* This file defines the user interface to the fann library.
   It is included from fixedfann.h, floatfann.h and doublefann.h and should
   NOT be included directly. If included directly it will react as if
   floatfann.h was included.
*/ 

/* Section: FANN Creation/Execution
   
   The FANN library is designed to be very easy to use. 
   A feedforward ann can be created by a simple <fann_create_standard> function, while
   other ANNs can be created just as easily. The ANNs can be trained by <fann_train_on_file>
   and executed by <fann_run>.
   
   All of this can be done without much knowledge of the internals of ANNs, although the ANNs created will
   still be powerfull and effective. If you have more knowledge about ANNs, and desire more control, almost
   every part of the ANNs can be parametized to create specialized and highly optimal ANNs.
 */
/* Group: Creation, Destruction & Execution */
	
#ifndef FANN_INCLUDE
/* just to allow for inclusion of fann.h in normal stuations where only floats are needed */ 
#ifdef FIXEDFANN
#include "fixedfann.h"
#else
#include "floatfann.h"
#endif	/* FIXEDFANN  */
	
#else
	
/* COMPAT_TIME REPLACEMENT */ 
#ifndef _WIN32
#include <sys/time.h>
#else	/* _WIN32 */
#if !defined(_MSC_EXTENSIONS) && !defined(_INC_WINDOWS)  
extern unsigned long __stdcall GetTickCount(void);

#else	/* _MSC_EXTENSIONS */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif	/* _MSC_EXTENSIONS */
#endif	/* _WIN32 */
		
#ifndef __fann_h__
#define __fann_h__
	
#ifdef __cplusplus
extern "C"
{
	
#ifndef __cplusplus
} /* to fool automatic indention engines */ 
#endif
#endif	/* __cplusplus */
 
#ifndef NULL
#define NULL 0
#endif	/* NULL */
 
/* ----- Macros used to define DLL external entrypoints ----- */ 
/*
 DLL Export, import and calling convention for Windows.
 Only defined for Microsoft VC++ FANN_EXTERNAL indicates
 that a function will be exported/imported from a dll
 FANN_API ensures that the DLL calling convention
 will be used for  a function regardless of the calling convention
 used when compiling.

 For a function to be exported from a DLL its prototype and
 declaration must be like this:
    FANN_EXTERNAL void FANN_API function(char *argument)

 The following ifdef block is a way of creating macros which
 make exporting from a DLL simple. All files within a DLL are
 compiled with the FANN_DLL_EXPORTS symbol defined on the
 command line. This symbol should not be defined on any project
 that uses this DLL. This way any other project whose source
 files include this file see FANN_EXTERNAL functions as being imported
 from a DLL, whereas a DLL sees symbols defined with this
 macro as being exported which makes calls more efficient.
 The __stdcall calling convention is used for functions in a
 windows DLL.

 The callback functions for fann_set_callback must be declared as FANN_API
 so the DLL and the application program both use the same
 calling convention.
*/ 
 
/*
 The following sets the default for MSVC++ 2003 or later to use
 the fann dll's. To use a lib or fixedfann.c, floatfann.c or doublefann.c
 with those compilers FANN_NO_DLL has to be defined before
 including the fann headers.
 The default for previous MSVC compilers such as VC++ 6 is not
 to use dll's. To use dll's FANN_USE_DLL has to be defined before
 including the fann headers.
*/ 
#if (_MSC_VER > 1300)
#ifndef FANN_NO_DLL
//#define FANN_USE_DLL
#endif	/* FANN_USE_LIB */
#endif	/* _MSC_VER */
#if defined(_MSC_VER) && (defined(FANN_USE_DLL) || defined(FANN_DLL_EXPORTS))
#ifdef FANN_DLL_EXPORTS
#define FANN_EXTERNAL __declspec(dllexport)
#else							/*  */
#define FANN_EXTERNAL __declspec(dllimport)
#endif	/* FANN_DLL_EXPORTS*/
#define FANN_API __stdcall
#else							/*  */
#define FANN_EXTERNAL
#define FANN_API
#endif	/* _MSC_VER */
/* ----- End of macros used to define DLL external entrypoints ----- */ 

#include "fann_error.h"
#include "fann_activation.h"
#include "fann_data.h"
#include "fann_internal.h"
#include "fann_train.h"
#include "fann_cascade.h"
#include "fann_io.h"

/* Function: fann_create_standard
	
	Creates a standard fully connected backpropagation neural network.

	There will be a bias neuron in each layer (except the output layer),
	and this bias neuron will be connected to all neurons in the next layer.
	When running the network, the bias nodes always emits 1.
	
	To destroy a <struct fann> use the <fann_destroy> function.

	Parameters:
		num_layers - The total number of layers including the input and the output layer.
		... - Integer values determining the number of neurons in each layer starting with the 
			input layer and ending with the output layer.
			
	Returns:
		A pointer to the newly created <struct fann>.
			
	Example:
		> // Creating an ANN with 2 input neurons, 1 output neuron, 
		> // and two hidden neurons with 8 and 9 neurons
		> struct fann *ann = fann_create_standard(4, 2, 8, 9, 1);
		
	See also:
		<fann_create_standard_array>, <fann_create_sparse>, <fann_create_shortcut>		
		
	This function appears in FANN >= 2.0.0.
*/ 
FANN_EXTERNAL struct fann *FANN_API fann_create_standard(unsigned int num_layers, ...);

/* Function: fann_create_standard_array
   Just like <fann_create_standard>, but with an array of layer sizes
   instead of individual parameters.

	Example:
		> // Creating an ANN with 2 input neurons, 1 output neuron, 
		> // and two hidden neurons with 8 and 9 neurons
		> unsigned int layers[4] = {2, 8, 9, 1};
		> struct fann *ann = fann_create_standard_array(4, layers);

	See also:
		<fann_create_standard>, <fann_create_sparse>, <fann_create_shortcut>

	This function appears in FANN >= 2.0.0.
*/ 
FANN_EXTERNAL struct fann *FANN_API fann_create_standard_array(unsigned int num_layers,
													           unsigned int *layers);

/* Function: fann_create_sparse

	Creates a standard backpropagation neural network, which is not fully connected.

	Parameters:
		connection_rate - The connection rate controls how many connections there will be in the
   			network. If the connection rate is set to 1, the network will be fully
   			connected, but if it is set to 0.5 only half of the connections will be set.
			A connection rate of 1 will yield the same result as <fann_create_standard>
		num_layers - The total number of layers including the input and the output layer.
		... - Integer values determining the number of neurons in each layer starting with the 
			input layer and ending with the output layer.
			
	Returns:
		A pointer to the newly created <struct fann>.

	See also:
		<fann_create_sparse_array>, <fann_create_standard>, <fann_create_shortcut>

	This function appears in FANN >= 2.0.0.
*/
FANN_EXTERNAL struct fann *FANN_API fann_create_sparse(float connection_rate, 
	                                                   unsigned int num_layers, ...);


/* Function: fann_create_sparse_array
   Just like <fann_create_sparse>, but with an array of layer sizes
   instead of individual parameters.

	See <fann_create_standard_array> for a description of the parameters.

	See also:
		<fann_create_sparse>, <fann_create_standard>, <fann_create_shortcut>

	This function appears in FANN >= 2.0.0.
*/
FANN_EXTERNAL struct fann *FANN_API fann_create_sparse_array(float connection_rate, 
	                                                         unsigned int num_layers, 
															 unsigned int *layers);

/* Function: fann_create_shortcut

	Creates a standard backpropagation neural network, which is not fully connected and which
	also has shortcut connections.

 	Shortcut connections are connections that skip layers. A fully connected network with shortcut 
	connections, is a network where all neurons are connected to all neurons in later layers. 
	Including direct connections from the input layer to the output layer.

	See <fann_create_standard> for a description of the parameters.

	See also:
		<fann_create_shortcut_array>, <fann_create_standard>, <fann_create_sparse>, 

	This function appears in FANN >= 2.0.0.
*/ 
FANN_EXTERNAL struct fann *FANN_API fann_create_shortcut(unsigned int num_layers, ...);

/* Function: fann_create_shortcut_array
   Just like <fann_create_shortcut>, but with an array of layer sizes
   instead of individual parameters.

	See <fann_create_standard_array> for a description of the parameters.

	See also:
		<fann_create_shortcut>, <fann_create_standard>, <fann_create_sparse>

	This function appears in FANN >= 2.0.0.
*/
FANN_EXTERNAL struct fann *FANN_API fann_create_shortcut_array(unsigned int num_layers,
															   unsigned int *layers);
/* Function: fann_destroy
   Destroys the entire network and properly freeing all the associated memmory.

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL void FANN_API fann_destroy(struct fann *ann);


/* Function: fann_run
	Will run input through the neural network, returning an array of outputs, the number of which being 
	equal to the number of neurons in the output layer.

	See also:
		<fann_test>

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL fann_type * FANN_API fann_run(struct fann *ann, fann_type * input);

/* Function: fann_randomize_weights
	Give each connection a random weight between *min_weight* and *max_weight*
   
	From the beginning the weights are random between -0.1 and 0.1.

	See also:
		<fann_init_weights>

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL void FANN_API fann_randomize_weights(struct fann *ann, fann_type min_weight,
												   fann_type max_weight);

/* Function: fann_init_weights
  	Initialize the weights using Widrow + Nguyen's algorithm.
	
 	This function behaves similarly to fann_randomize_weights. It will use the algorithm developed 
	by Derrick Nguyen and Bernard Widrow to set the weights in such a way 
	as to speed up training. This technique is not always successful, and in some cases can be less 
	efficient than a purely random initialization.

	The algorithm requires access to the range of the input data (ie, largest and smallest input), 
	and therefore accepts a second argument, data, which is the training data that will be used to 
	train the network.

	See also:
		<fann_randomize_weights>, <fann_read_train_from_file>

	This function appears in FANN >= 1.1.0.
*/ 
FANN_EXTERNAL void FANN_API fann_init_weights(struct fann *ann, struct fann_train_data *train_data);

/* Function: fann_print_connections
	Will print the connections of the ann in a compact matrix, for easy viewing of the internals 
	of the ann.

	The output from fann_print_connections on a small (2 2 1) network trained on the xor problem
	>Layer / Neuron 012345
	>L   1 / N    3 BBa...
	>L   1 / N    4 BBA...
	>L   1 / N    5 ......
	>L   2 / N    6 ...BBA
	>L   2 / N    7 ......
		  
	This network have five real neurons and two bias neurons. This gives a total of seven neurons 
	named from 0 to 6. The connections between these neurons can be seen in the matrix. "." is a 
	place where there is no connection, while a character tells how strong the connection is on a 
	scale from a-z. The two real neurons in the hidden layer (neuron 3 and 4 in layer 1) has 
	connection from the three neurons in the previous layer as is visible in the first two lines. 
	The output neuron (6) has connections form the three neurons in the hidden layer 3 - 5 as is 
	visible in the fourth line.

	To simplify the matrix output neurons is not visible as neurons that connections can come from, 
	and input and bias neurons are not visible as neurons that connections can go to.

	This function appears in FANN >= 1.2.0.
*/ 
FANN_EXTERNAL void FANN_API fann_print_connections(struct fann *ann);

/* Group: Parameters */
/* Function: fann_print_parameters

  	Prints all of the parameters and options of the ANN 

	This function appears in FANN >= 1.2.0.
*/ 
FANN_EXTERNAL void FANN_API fann_print_parameters(struct fann *ann);


/* Function: fann_get_num_input

   Get the number of input neurons.

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL unsigned int FANN_API fann_get_num_input(struct fann *ann);


/* Function: fann_get_num_output

   Get the number of output neurons.

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL unsigned int FANN_API fann_get_num_output(struct fann *ann);


/* Function: fann_get_total_neurons

   Get the total number of neurons in the entire network. This number does also include the 
	bias neurons, so a 2-4-2 network has 2+4+2 +2(bias) = 10 neurons.

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL unsigned int FANN_API fann_get_total_neurons(struct fann *ann);


/* Function: fann_get_total_connections

   Get the total number of connections in the entire network.

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL unsigned int FANN_API fann_get_total_connections(struct fann *ann);

#ifdef FIXEDFANN
	
/* Function: fann_get_decimal_point

	Returns the position of the decimal point in the ann.

	This function is only available when the ANN is in fixed point mode.

	The decimal point is described in greater detail in the tutorial <Fixed Point Usage>.

	See also:
		<Fixed Point Usage>, <fann_get_multiplier>, <fann_save_to_fixed>, <fann_save_train_to_fixed>

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL unsigned int FANN_API fann_get_decimal_point(struct fann *ann);


/* Function: fann_get_multiplier

    returns the multiplier that fix point data is multiplied with.

	This function is only available when the ANN is in fixed point mode.

	The multiplier is the used to convert between floating point and fixed point notation. 
	A floating point number is multiplied with the multiplier in order to get the fixed point
	number and visa versa.

	The multiplier is described in greater detail in the tutorial <Fixed Point Usage>.

	See also:
		<Fixed Point Usage>, <fann_get_decimal_point>, <fann_save_to_fixed>, <fann_save_train_to_fixed>

	This function appears in FANN >= 1.0.0.
*/ 
FANN_EXTERNAL unsigned int FANN_API fann_get_multiplier(struct fann *ann);

#endif	/* FIXEDFANN */

#ifdef __cplusplus
#ifndef __cplusplus
/* to fool automatic indention engines */ 
{
	
#endif
} 
#endif	/* __cplusplus */
	
#endif	/* __fann_h__ */
	
#endif /* NOT FANN_INCLUDE */
