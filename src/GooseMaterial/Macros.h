/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

================================================================================================= */

#ifndef GOOSEMATERIAL_MACROS_H
#define GOOSEMATERIAL_MACROS_H

#define GOOSEMATERIAL_WORLD_VERSION 0
#define GOOSEMATERIAL_MAJOR_VERSION 0
#define GOOSEMATERIAL_MINOR_VERSION 7

#define GOOSEMATERIAL_VERSION_AT_LEAST(x,y,z) \
  (GOOSEMATERIAL_WORLD_VERSION>x || (GOOSEMATERIAL_WORLD_VERSION>=x && \
  (GOOSEMATERIAL_MAJOR_VERSION>y || (GOOSEMATERIAL_MAJOR_VERSION>=y && \
                                     GOOSEMATERIAL_MINOR_VERSION>=z))))

#define GOOSEMATERIAL_VERSION(x,y,z) \
  (GOOSEMATERIAL_WORLD_VERSION==x && \
   GOOSEMATERIAL_MAJOR_VERSION==y && \
   GOOSEMATERIAL_MINOR_VERSION==z)

#endif
