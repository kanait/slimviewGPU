////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Mesh class for Rendering
//
// Copyright (c) 2003-2010 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _MESHR_HXX
#define _MESHR_HXX 1

#include "envDep.h"
#include "mydef.h"

#if defined(_WINDOWS)
#include "stdafx.h"
#endif

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
using namespace std;

#include <Point3.h>
#include <Vector3.h>
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

// #include "MeshL.hxx"

#define ASSIGN_VERTEX 0
#define ASSIGN_FACE   1

class MeshR {

public:
  
  MeshR() { clear(); };
  ~MeshR() { clear(); };

  void clear() {

    points_.clear();
    normals_.clear();
    texcoords_.clear();
    colors_.clear();
    indices_.clear();
    nindices_.clear();
    fnormals_.clear();
    face_mates_.clear();
    s_begin_.clear();
    s_index_size_.clear();

    n_points_ = 0;
    n_normals_ = 0;
    n_texcoords_ = 0;
    n_colors_ = 0;
    n_indices_ = 0;
    n_nindices_ = 0;
    n_tex_ = 2; // x, y coordinates

    isNormalized_ = false;
  };

//   void setMesh( MeshL& );

  bool empty() const { 
    if ( numPoints() || numNormals() || numTexcoords() || numColors() || numFaces() )
      return false;
    return true;
  };

  // get real number of elements
  unsigned int numPoints() const { return (unsigned int) (points_size() / 3.0f); };
  unsigned int numNormals() const { return (unsigned int) (normals_size() / 3.0f); };
  unsigned int numTexcoords() const { return (unsigned int) ((float) texcoords_size() / (float) n_tex_); };
  unsigned int numColors() const { return (unsigned int) (colors_size() / 3.0f); };
  unsigned int numFaces() const { return (unsigned int) (indices_size() / 3.0f); };

  // get the number of internal elements
  unsigned int points_size() const { return n_points_; };
  unsigned int normals_size() const { return n_normals_; };
  unsigned int texcoords_size() const { return n_texcoords_; };
  unsigned int colors_size() const { return n_colors_; };
  unsigned int indices_size() const { return n_indices_; };
  unsigned int nindices_size() const { return n_nindices_; };

  std::vector<float>& points() { return points_; };
  std::vector<float>& normals() { return normals_; };
  std::vector<float>& texcoords() { return texcoords_; };
  std::vector<unsigned char>& colors() { return colors_; };
  std::vector<unsigned int>& color_ids() { return color_ids_; };
  std::vector<unsigned int>& indices() { return indices_; };
  std::vector<unsigned int>& nindices() { return nindices_; };
  std::vector<float>& fnormals() { return fnormals_; };
  std::vector<int>& face_mates() { return face_mates_; };
  std::vector<unsigned int>& s_begin() { return s_begin_; };
  std::vector<unsigned int>& s_index_size() { return s_index_size_; };

//  Point3f& point( int i ) { return Point3f( points_[3*i], points_[3*i+1], points_[3*i+2] ); };
//  Vector3f& normal( int i ) { return Vector3f( normals_[3*i], normals_[3*i+1], normals_[3*i+2] ); };
//  Point3f& texcoord( int i ) { return Point3f( texcoords_[3*i], texcoords_[3*i+1], texcoords_[3*i+2] ); };
//  Point3f& color( int i ) { return Point3f( colors_[3*i], colors_[3*i+1], colors_[3*i+2] ); };
  int index( int i ) const { return indices_[i]; };
  int nindex( int i ) const { return nindices_[i]; };
  float fnormal( int i ) const { return fnormals_[i]; };
  int face_mate( int i ) const { return face_mates_[i]; };

  int n_tex() const { return n_tex_; };
  void setNTex( int n ) { n_tex_ = n; };

  void setColorAssigned( unsigned short n ) { color_assigned_ = n; };
  unsigned short colorAssigned() const { return color_assigned_; };

  //
  // "add" elements
  //

  void addIndex( unsigned int i ) {
    indices_.push_back( i ); n_indices_++;
  };

  void addNIndex( unsigned int i ) {
    nindices_.push_back( i ); n_nindices_++;
  };

  void addPoint( float x, float y, float z ) {
    points_.push_back( x ); n_points_++;
    points_.push_back( y ); n_points_++;
    points_.push_back( z ); n_points_++;
  };

  void addNormal( float x, float y, float z ) {
    normals_.push_back( x ); n_normals_++;
    normals_.push_back( y ); n_normals_++;
    normals_.push_back( z ); n_normals_++;
  };

  void addTexcoord( float x, float y ) {
    n_tex_ = 2;
    texcoords_.push_back( x ); n_texcoords_++;
    texcoords_.push_back( y ); n_texcoords_++;
  };

  void addTexcoord( float x, float y, float z ) {
    n_tex_ = 3;
    texcoords_.push_back( x ); n_texcoords_++;
    texcoords_.push_back( y ); n_texcoords_++;
    texcoords_.push_back( z ); n_texcoords_++;
  };

  void addColor( unsigned char x, unsigned char y, unsigned char z ) {
    colors_.push_back( x ); n_colors_++;
    colors_.push_back( y ); n_colors_++;
    colors_.push_back( z ); n_colors_++;
  };
  
  //
  // "set" elements
  //

  void setPoint( int i, Point3f& p ) {
    points_[i]   = p.x;
    points_[i+1] = p.y;
    points_[i+2] = p.z;
  };

  void setPoint( int i, float x, float y, float z ) {
    points_[i]   = x;
    points_[i+1] = y;
    points_[i+2] = z;
  };

  void setNormal( int i, Vector3f& p ) {
    normals_[i]   = p.x;
    normals_[i+1] = p.y;
    normals_[i+2] = p.z;
  };

  void setNormal( int i, float x, float y, float z ) {
    normals_[i]   = x;
    normals_[i+1] = y;
    normals_[i+2] = z;
  };
    
  void setTexcoord( int i, Point3f& p ) {
    texcoords_[i]   = p.x;
    texcoords_[i+1] = p.y;
    texcoords_[i+2] = p.z;
  };

  void setTexcoord( int i, float x, float y, float z ) {
    texcoords_[i]   = x;
    texcoords_[i+1] = y;
    texcoords_[i+2] = z;
  };

  void setTexcoord( int i, float x, float y ) {
    texcoords_[i]   = x;
    texcoords_[i+1] = y;
  };

//   void setColor( int i, Point3& p ) {
//     colors_[i]   = p.x;
//     colors_[i+1] = p.y;
//     colors_[i+2] = p.z;
//   };

  void setColor( int i, unsigned char x, unsigned char y, unsigned char z ) {
    colors_[i]   = x;
    colors_[i+1] = y;
    colors_[i+2] = z;
  };

  void setColorId( int i, unsigned int n ) {
    color_ids_[i] = n;
  };

  void setIndex( int i, unsigned int f ) {
    indices_[i] = f;
  };

  void setNIndex( int i, unsigned int f ) {
    nindices_[i] = f;
  };

  void deleteTexcoords() { texcoords_.clear(); n_texcoords_ = 0; };

  // reserve
  void reservePoints( int n )  { n_points_  = n * nXYZ;     points_.resize( n_points_ ); };
  void reserveNormals( int n ) { n_normals_ = n * nXYZ;     normals_.resize( n_normals_ ); };
  void reserveTexcoords( int n, int t ) { n_texcoords_ = n * t; texcoords_.resize( n_texcoords_ ); };
  void reserveColors( int n ) { n_colors_ = n * nXYZ;       colors_.resize( n_colors_ ); };
  void reserveColorIds( int n ) { n_color_ids_ = n;         color_ids_.resize( n_color_ids_ ); };
  void reserveIndices( int n ) { n_indices_ = n * TRIANGLE; indices_.resize( n_indices_ ); };
  void reserveNIndices( int n ) { n_nindices_ = n * TRIANGLE; nindices_.resize( n_nindices_ ); };

  // resize
  void resizePoints( int n )  { n_points_  += (n * nXYZ);     points_.resize( n_points_ ); };
  void resizeNormals( int n ) { n_normals_ += (n * nXYZ);     normals_.resize( n_normals_ ); };
  void resizeIndices( int n ) { n_indices_ += (n * TRIANGLE); indices_.resize( n_indices_ ); };
  
  void getFacePoints( unsigned int face_id, Point3f& p0, Point3f& p1, Point3f& p2 ) {
    unsigned int i0 = nXYZ * indices_[ TRIANGLE * face_id ];
    unsigned int i1 = nXYZ * indices_[ TRIANGLE * face_id + 1 ];
    unsigned int i2 = nXYZ * indices_[ TRIANGLE * face_id + 2 ];
    p0.set( points_[ i0 ], points_[ i0 + 1 ], points_[ i0 + 2 ] );
    p1.set( points_[ i1 ], points_[ i1 + 1 ], points_[ i1 + 2 ] );
    p2.set( points_[ i2 ], points_[ i2 + 1 ], points_[ i2 + 2 ] );
  };

  void getFaceMates( unsigned int face_id, unsigned int m[TRIANGLE] ) {
    m[0] = (int) ( (float) face_mates_[ TRIANGLE * face_id ] / 3.0f );
    m[1] = (int) ( (float) face_mates_[ TRIANGLE * face_id + 1 ] / 3.0f );
    m[2] = (int) ( (float) face_mates_[ TRIANGLE * face_id + 2 ] / 3.0f );
  };

  void calcNormal( Point3f& p0, Point3f& p1, Point3f& p2, Vector3f& nrm ) {
    Vector3f v1( p1 - p0 );
    Vector3f v2( p2 - p0 );
    nrm.cross(v1, v2);
    nrm.normalize();
  };
  
  float calcFaceArea( Point3f& p0, Point3f& p1, Point3f& p2 ) {
    Vector3f v1( p1 - p0 );
    Vector3f v2( p2 - p0 );
    Vector3f a;
    a.cross(v1, v2);
    return a.length() * .5f;
  };
  
  void scale( float f ) {
    for ( int i = 0; i < points_.size(); i += nXYZ )
      {
	points_[i] *= f;
	points_[i+1] *= f;
	points_[i+2] *= f;
      }
  };


  bool isNormalized() const { return isNormalized_; };
  void setIsNormalized( bool t ) { isNormalized_ = t; };

  void normalize() {
    if ( isNormalized_ ) return;

    std::cout << "normalize ... " << endl;

    getTransScale();

    for ( int i = 0; i < points_.size(); i += nXYZ )
      {
	float x = (points_[i]   - center_.x) / maxlen_;
	float y = (points_[i+1] - center_.y) / maxlen_;
	float z = (points_[i+2] - center_.z) / maxlen_;
	points_[i]   = x;
	points_[i+1] = y;
	points_[i+2] = z;
      }

    isNormalized_ = true;

    std::cout << "done." << endl;
  };

  void normalize( Point3f& center, float maxlen ) {
    if ( isNormalized_ ) return;

    std::cout << "normalize ... " << endl;

    for ( int i = 0; i < points_.size(); i += nXYZ )
      {
	float x = (points_[i]   - center.x) / maxlen;
	float y = (points_[i+1] - center.y) / maxlen;
	float z = (points_[i+2] - center.z) / maxlen;
	points_[i]   = x;
	points_[i+1] = y;
	points_[i+2] = z;
      }

    isNormalized_ = true;

    std::cout << "done." << endl;
  };

  void getTransScale() {
    getTransScale( vmin_, vmax_, center_, length_, &maxlen_ );
  };
  void getTransScale( Point3f& vmin, Point3f& vmax, Point3f& cen,
		      Point3f& len, float* maxlen ) {
    for ( int i = 0; i < points_.size(); i += nXYZ )
      {
	float x = points_[i]; 
	float y = points_[i+1]; 
	float z = points_[i+2];
	if ( i )
	  {
	    if (x > vmax.x) vmax.x = x;
	    if (x < vmin.x) vmin.x = x;
	    if (y > vmax.y) vmax.y = y;
	    if (y < vmin.y) vmin.y = y;
	    if (z > vmax.z) vmax.z = z;
	    if (z < vmin.z) vmin.z = z;
	  }
	else
	  {
	    vmax.set( x, y, z ); vmin.set( x, y, z );
	  }
      }

    cen = vmax + vmin; cen.scale(.5);
    len = vmax - vmin;
    *maxlen = (std::fabs(len.x) > std::fabs(len.y) )
      ? std::fabs(len.x) : std::fabs(len.y); 
    *maxlen = ( *maxlen > std::fabs(len.z) ) ? *maxlen : std::fabs(len.z);
  };

  Point3f& center() { return center_; };
  Point3f& vmin() { return vmin_; };
  Point3f& vmax() { return vmax_; };
  Point3f& length() { return length_; };
  float maxlen() const { return maxlen_; }; 
  
  void createFaceNormals() {
    if ( !fnormals_.empty() ) return;
  
    fnormals_.resize( numFaces() * nXYZ );

    for ( int i = 0; i < numFaces(); ++i )
      {
	Vector3f nrm;
	calcNormal( Point3f( points_[ nXYZ * indices_[ TRIANGLE * i ] ],
			     points_[ nXYZ * indices_[ TRIANGLE * i ] + 1 ],
			     points_[ nXYZ * indices_[ TRIANGLE * i ] + 2 ] ),
		    Point3f( points_[ nXYZ * indices_[ TRIANGLE * i + 1 ] ],
			     points_[ nXYZ * indices_[ TRIANGLE * i + 1 ] + 1 ],
			     points_[ nXYZ * indices_[ TRIANGLE * i + 1 ] + 2 ] ),
		    Point3f( points_[ nXYZ * indices_[ TRIANGLE * i + 2 ] ],
			     points_[ nXYZ * indices_[ TRIANGLE * i + 2 ] + 1 ],
			     points_[ nXYZ * indices_[ TRIANGLE * i + 2 ] + 2 ] ),
		    nrm );
	fnormals_[ nXYZ * i ] = nrm.x;
	fnormals_[ nXYZ * i + 1 ] = nrm.y;
	fnormals_[ nXYZ * i + 2 ] = nrm.z;
      }
  };

  void createVertexNormals() {
    std::cout << "calculate normals ... " << std::endl;

    //   std::cout << "points " << points_.size() << std::endl;

    normals_.resize( points_.size() );
  
    for ( int i = 0; i < normals_.size(); ++i ) normals_[i] = .0f;
    std::vector<float> weights; weights.resize(  numPoints() );
    for ( int i = 0; i < weights.size(); ++i ) weights[i] = .0f;

    // store face normals multiplied by area weights
    for ( int i = 0; i < indices_.size(); i += TRIANGLE )
      {
	unsigned int id0 = indices_[i];
	unsigned int id1 = indices_[i+1];
	unsigned int id2 = indices_[i+2];

	Point3f p0 ( points_[ nXYZ * id0 ], 
		     points_[ nXYZ * id0 + 1 ],
		     points_[ nXYZ * id0 + 2 ] );
	Point3f p1 ( points_[ nXYZ * id1 ], 
		     points_[ nXYZ * id1 + 1 ],
		     points_[ nXYZ * id1 + 2 ] );
	Point3f p2 ( points_[ nXYZ * id2 ], 
		     points_[ nXYZ * id2 + 1 ],
		     points_[ nXYZ * id2 + 2 ] );
      
	Vector3f nrm;
	calcNormal( p0, p1, p2, nrm );
	float w = calcFaceArea( p0, p1, p2 );
	nrm.scale( w );

	normals_[ nXYZ * (id0) ]     += nrm.x;
	normals_[ nXYZ * (id0) + 1 ] += nrm.y;
	normals_[ nXYZ * (id0) + 2 ] += nrm.z;
	weights[ id0 ] += w;
      
	normals_[ nXYZ * (id1) ]     += nrm.x;
	normals_[ nXYZ * (id1) + 1 ] += nrm.y;
	normals_[ nXYZ * (id1) + 2 ] += nrm.z;
	weights[ id1 ] += w;

	normals_[ nXYZ * (id2) ]     += nrm.x;
	normals_[ nXYZ * (id2) + 1 ] += nrm.y;
	normals_[ nXYZ * (id2) + 2 ] += nrm.z;
	weights[ id2 ] += w;
      }

    // divided by weight and normalize
    for ( int i = 0; i < numPoints(); ++i )
      {
	normals_[ nXYZ * i ]     /= weights[i];
	normals_[ nXYZ * i + 1 ] /= weights[i];
	normals_[ nXYZ * i + 2 ] /= weights[i];

	Vector3f nrm( normals_[ nXYZ * i ], normals_[ nXYZ * i + 1 ], normals_[ nXYZ * i + 2 ] );
	nrm.normalize();

	normals_[ nXYZ * i ]     = nrm.x;
	normals_[ nXYZ * i + 1 ] = nrm.y;
	normals_[ nXYZ * i + 2 ] = nrm.z;
      }

    n_normals_ = normals_.size();

    std::cout << "done. n " << (int) (n_normals_ / 3.0) << std::endl;
  };

  void createVertexNormalsWithSF() {
    if ( fnormals_.empty() ) createFaceNormals();
  
    unsigned int id_vn = 0;
    // ���_���̖@���O���[�v�̐�
    std::vector< std::vector<int> > n_vn( numPoints() );
    // ���_���̊e�@���O���[�v�̒��_ID
    std::vector< std::vector< std::vector<unsigned int> > > vid_vn( numPoints() );
    // ���_���̊e�@���O���[�v�̖@���i�O���[�v���ŕ��ω����Ă���j
    std::vector< std::vector<Vector3f> > vnrm( numPoints() );

    // ���_���Ƃɖ@����ǉ�
    // ����臒l���傫���p�x�̖@���͕ʂɊǗ�
    for ( int i = 0; i < numFaces(); ++i )
      {
	unsigned int id[3];
	id[0] = indices_[TRIANGLE * i];
	id[1] = indices_[TRIANGLE * i + 1];
	id[2] = indices_[TRIANGLE * i + 2];

	// ���ʂ̖@��
	Vector3f nf( fnormals_[ nXYZ * i ], 
		     fnormals_[ nXYZ * i + 1 ], 
		     fnormals_[ nXYZ * i + 2 ] );

	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    bool flag = false;
	    // 	  cout << "j = " << j << endl;
	    // ���_�̖@���O���[�v�̐�
	    for ( int k = 0; k < n_vn[id[j]].size(); ++k )
	      {
		// 	      cout << id[j] << " " << n_vn[id[j]][k] << vnrm[id[j]][k] << endl;
		// �i�[����Ă���O���[�v�̖@��
		Vector3f& nf0 = vnrm[id[j]][k];
#if 0
		cout << "id[j] " << id[j] << " [k] " << k << " vn " << " size " << n_vn[id[j]][k] << endl;
		for ( int l = 0; l < n_vn[id[j]][k]; ++l )
		  {
		    cout << "\t vid_vn " << vid_vn[id[j]][k][l] << endl;
		  }
#endif
		// 臒l�ȓ���������@���𕽋ς���
		if ( nf.angle( nf0 ) < (M_PI / 9.0f) ) 
		  {
		    unsigned int kvn = n_vn[id[j]][k];
		    nf0.scale( (float) kvn );
		    nf0.add( nf );
		    nf0.scale( 1.0f / (float) (kvn+1) );
		    ++(n_vn[id[j]][k]);
		    vid_vn[id[j]][k].push_back( id_vn++ );

		    flag = true;
		    break;
		  }
	      }

	    // �V���ɖ@���O���[�v��ǉ�
	    if ( flag == false )
	      {
		// 	      cout << "aa" << endl;
		n_vn[id[j]].push_back(1);
		vnrm[id[j]].push_back( nf );
		//vid_vn[id[j]].push_back( 1 );
		std::vector<unsigned int> nvn(1); nvn[0] = id_vn++;
		vid_vn[id[j]].push_back( nvn );
	      }
	  }
      }

    // �@���x�N�g���̊i�[
    if ( !normals_.empty() ) normals_.clear();
    normals_.resize( numFaces() * TRIANGLE * nXYZ );

    for ( int i = 0; i < numFaces(); ++i )
      {
	unsigned int id[3];
	id[0] = indices_[TRIANGLE * i];
	id[1] = indices_[TRIANGLE * i + 1];
	id[2] = indices_[TRIANGLE * i + 2];

	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    for ( int k = 0; k < n_vn[id[j]].size(); ++k )
	      {
		Vector3f& nf = vnrm[id[j]][k];
		nf.normalize();
		for ( int l = 0; l < n_vn[id[j]][k]; ++l )
		  {
		    unsigned int nid = vid_vn[id[j]][k][l];
		    // 		  cout << "nid " << nid << endl;
		    normals_[ nXYZ * nid ] = nf.x;
		    normals_[ nXYZ * nid + 1 ] = nf.y;
		    normals_[ nXYZ * nid + 2 ] = nf.z;
		  }
		// 	      cout << endl;
	      }
	  }

      }
  };

  // �ʂ̃��C�g�Q�̐���
  void createFaceMates() {
    std::cout << "create index pair map ... " << std::endl;

    typedef std::pair< unsigned int, unsigned int > type_uu;
    typedef std::pair< unsigned int, std::pair< unsigned int, unsigned int> > type_up;
    typedef std::multimap< unsigned int, std::pair< unsigned int, unsigned int> >::iterator type_upIter;
    // multimap< id0, pair< id1, fid > >
    std::multimap< unsigned int, std::pair< unsigned int, unsigned int> > indices_pair;

    // �̈�̊m��
    face_mates_.resize( indices_.size() );

    for( int i = 0; i < indices_.size(); i += TRIANGLE )
      {
#if 1
	unsigned int id[3];
	id[0] = indices_[ i ];
	id[1] = indices_[ i + 1 ];
	id[2] = indices_[ i + 2 ];

	indices_pair.insert( type_up(id[0], type_uu(id[1], i)) );
	indices_pair.insert( type_up(id[1], type_uu(id[2], i)) );
	indices_pair.insert( type_up(id[2], type_uu(id[0], i)) );
#endif

#if 0
	unsigned int id[3];
	id[0] = indices_[ i ];
	id[1] = indices_[ i + 1 ];
	id[2] = indices_[ i + 2 ];

	type_upIter cIterI[3];
	cIterI[0] = indices_pair.insert( type_up(id[0], type_uu(id[1], i)) );
	cIterI[1] = indices_pair.insert( type_up(id[1], type_uu(id[2], i)) );
	cIterI[2] = indices_pair.insert( type_up(id[2], type_uu(id[0], i)) );

	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    unsigned int id0 = indices_[ i + j ];
	    unsigned int id1 = ( j != 2 ) ? indices_[ i + (j + 1) ] : indices_[ i ];

	    // id1 ���L�[�Ƃ��� pair ��������
	    std::pair< type_upIter, type_upIter > cIterPair
	      = indices_pair.equal_range( id1 );
	  
	    face_mates_[ i + j ] = -1;
	    bool found = false;
	    for ( type_upIter cIter = cIterPair.first; cIter != cIterPair.second; ++cIter )
	      {
		// pair �̂��� �ŏ��̃L�[�� id0 �̂��̂�������
		if ( (*cIter).second.first == id0 ) // found!
		  {
		    found = true;
		    face_mates_[ i + j ] = (*cIter).second.second;
		    face_mates_[ (*cIter).second.second ] = i;

		    indices_pair.erase( cIter );
		    indices_pair.erase( cIterI[j] );
		    break;
		  }
	      }
#if 0
	    if ( found == false )
	      {
		indices_pair.insert( type_up(id0, type_uu(id1, i)) );
	      }
#endif
	  }
#endif
      }

#if 1
    // ���C�g�̍쐬
    std::cout << "create face mates ... " << std::endl;
  
    for( int i = 0; i < indices_.size(); i += TRIANGLE )
      {
	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    unsigned int id0 = indices_[ i + j ];
	    unsigned int id1 = ( j != 2 ) ? indices_[ i + (j + 1) ] : indices_[ i ];

	    // id1 ���L�[�Ƃ��� pair ��������
	    std::pair< type_upIter, type_upIter > cIterPair
	      = indices_pair.equal_range( id1 );
	  
	    face_mates_[ i + j ] = -1;
	    for ( type_upIter cIter = cIterPair.first; cIter != cIterPair.second; ++cIter )
	      {
		// pair �̂��� �ŏ��̃L�[�� id0 �̂��̂�������
		if ( (*cIter).second.first == id0 ) // found!
		  {
		    face_mates_[ i + j ] = (*cIter).second.second;
		  }
	      }
	  }
      }

#endif

    std::cout << "done. " << std::endl;
  };

  void deleteFaceMates() { face_mates_.clear(); };

  void printInfo() {
    std::cout << "mesh " << " ";
    if ( s_begin_.size() ) std::cout << " s " << s_begin_.size() << " ";
    if ( points_.size() ) std::cout << " v " << numPoints() << " ";
    if ( normals_.size() ) std::cout << " n " << numNormals() << " ";
    if ( colors_.size() ) std::cout << " c " << numColors() << " ";
    if ( texcoords_.size() ) std::cout << " t " << numTexcoords() << " ";
    if ( indices_.size() ) std::cout << " f " << numFaces() << " ";
    std::cout << std::endl;
  };

  void getBBox( Point3f& vmax, Point3f& vmin ) {

    for ( int i = 0; i < points_.size(); i += 3 )
      {
	Point3f p( points_[i], points_[i+1], points_[i+2] );
	if ( i )
	  {
	    if (p.x > vmax.x) vmax.x = p.x;
	    if (p.x < vmin.x) vmin.x = p.x;
	    if (p.y > vmax.y) vmax.y = p.y;
	    if (p.y < vmin.y) vmin.y = p.y;
	    if (p.z > vmax.z) vmax.z = p.z;
	    if (p.z < vmin.z) vmin.z = p.z;
	  }
	else
	  {
	    vmax.set( p ); vmin.set( p );
	  }
      }
  };
  
private:

  // vertex points
  unsigned int n_points_;
  std::vector<float> points_;

  // normals
  unsigned int n_normals_;
  std::vector<float> normals_;
 
  // texture coordinates
  unsigned int n_tex_;
  unsigned int n_texcoords_;
  std::vector<float> texcoords_;

  // colors
  unsigned short color_assigned_;
  unsigned int n_colors_;
  std::vector<unsigned char> colors_;

  // assigned color for vertices or faces
  unsigned int n_color_ids_;
  std::vector<unsigned int> color_ids_;

  // face indices
  unsigned int n_indices_;
  std::vector<unsigned int> indices_;

  // face normal indices
  unsigned int n_nindices_;
  std::vector<unsigned int> nindices_;

//   int n_indices_;
//   std::vector<unsigned int> normal_indices_;

  // face normals
  std::vector<float> fnormals_;

  // face mates
  std::vector<int> face_mates_;

  // solid
  std::vector<unsigned int> s_begin_;
  std::vector<unsigned int> s_index_size_;

  // for object
  Point3f vmin_;
  Point3f vmax_;
  Point3f center_;
  Point3f length_;
  float   maxlen_;

  // normalized
  bool isNormalized_;

};

#endif // _MESHR_HXX