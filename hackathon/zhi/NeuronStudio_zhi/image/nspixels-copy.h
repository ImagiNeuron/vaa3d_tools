#ifndef __NS_IMAGE_PIXELS_COPY_H__
#define __NS_IMAGE_PIXELS_COPY_H__

#include <image/nspixels.h>

NS_DECLS_BEGIN

enum
	{
	NS_PIXEL_PROC_COPY_PIXEL_TYPE,
	NS_PIXEL_PROC_COPY_SRC_PIXELS,
	NS_PIXEL_PROC_COPY_WIDTH,
	NS_PIXEL_PROC_COPY_HEIGHT,
	NS_PIXEL_PROC_COPY_LENGTH,
	NS_PIXEL_PROC_COPY_ROW_ALIGN,
	NS_PIXEL_PROC_COPY_DEST_PIXELS,
	NS_PIXEL_PROC_COPY_REGION,
	NS_PIXEL_PROC_COPY_PROGRESS,

	NS_PIXEL_PROC_COPY_NUM_PARAMS
	};

NS_IMPEXP NsProc* ns_pixel_proc_copy( void );

NS_DECLS_END

#endif/* __NS_IMAGE_PIXELS_COPY_H__ */
