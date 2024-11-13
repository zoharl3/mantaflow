/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * QT OpenGL widget
 *
 ******************************************************************************/

#include "glwidget.h"
#ifdef __APPLE__
#   include <OpenGL/glu.h>
#else
#   include <GL/glu.h>
#endif
#include <cmath>
#include "painter.h"

namespace Manta {

GLWidget::GLWidget(QWidget* p): QGLWidget(QGLFormat(QGL::SampleBuffers), p), mRotX(0), mRotY(0), mGridsize(0), mScreenshotNumber(0),
		mWidth(800), mHeight(600)
{
    m_b2D = false;
	mPlaneDim = 2; // Y plane
	mPlane = -1;
	mCamPos = Vec3(0, 0, -1.5);
	for (int i=0; i<MoveDirNum; i++) 
		mMoveState[i] = false;
	mMoveFast = false;
	
	setAutoBufferSwap(true);
	setFocusPolicy(Qt::ClickFocus);
	startTimer(10);
}

GLWidget::~GLWidget()
{

}

QSize GLWidget::minimumSizeHint() const
{
	return QSize(400, 300);
}

QSize GLWidget::sizeHint() const
{
	return QSize(mWidth, mHeight);
}
void GLWidget::windowSize(int w, int h) {
	mWidth = w; mHeight = h;
}

void GLWidget::initializeGL()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();     
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearDepth(1.0);    
}

void GLWidget::paintGL()
{
	if (mGridsize.max() == 0) return;
	glDepthFunc(GL_ALWAYS);
	float c = 0; // 0(default), 1(fig)
    glClearColor( c, c, c, 1.0f ); // zl background color
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable( GL_DEPTH_TEST );
	//glEnable(GL_POLYGON_OFFSET_FILL);
	//glPolygonOffset(0,0);
	
	// camera transform
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef( mCamPos.x, mCamPos.y, mCamPos.z );
	glRotatef( mRotX, 1., 0., 0. );
	glRotatef(mRotY,  0.,1.,0.);
	Real dx = 1.0 / (Real)mGridsize.max();
	Vec3 sz = toVec3(mGridsize) * (-0.5f * dx);
	
	glTranslatef( sz.x, sz.y, sz.z );
	Q_EMIT paintSub();
}

void GLWidget::resizeGL(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if ( m_b2D ) {
        GLfloat fov = 45;
        GLfloat zNear = 1.3f;
        GLfloat zFar = 100.0f;
        GLfloat aspect = float( w ) / float( h );
        GLfloat fH = tan( float( fov / 360.0f * M_PI ) ) * zNear;
        GLfloat fW = fH * aspect;
        glOrtho( -fW, fW, -fH, fH, zNear, zFar );
    } else {
        GLfloat fov = 45;
        GLfloat zNear = 0.05f;
        GLfloat zFar = 100.0f;
        GLfloat aspect = float( w ) / float( h );
        GLfloat fH = tan( float( fov / 360.0f * M_PI ) ) * zNear;
        GLfloat fW = fH * aspect;
        glFrustum( -fW, fW, -fH, fH, zNear, zFar );
    }

	glMatrixMode(GL_MODELVIEW);
}

void GLWidget::mouseReleaseEvent(QMouseEvent* event) {
	// only do tooltip if not moving
	QPoint pos = event->pos();
	if ((mDownPos - pos).manhattanLength() == 0) {
		// get GL transform matrices
		int viewport[4];
		GLdouble modelMatrix[16], projMatrix[16];
		glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
		glGetIntegerv(GL_VIEWPORT,viewport);
		
		// obtain click line
		GLdouble line[6], wx=pos.x(), wy=viewport[3]-pos.y();
		if (!gluUnProject(wx,wy,0,modelMatrix,projMatrix,viewport,&line[0],&line[1],&line[2])) return;
		if (!gluUnProject(wx,wy,1.0,modelMatrix,projMatrix,viewport,&line[3],&line[4],&line[5])) return;
		
		// calculate intersection with plane
		Q_EMIT clickLine(event->globalPos(), line[0],line[1],line[2],line[3],line[4],line[5]);
	}
}

void GLWidget::mouseMoveEvent(QMouseEvent* e)
{
	const float speedRot = 0.2f, speedPan = 0.002f;
	
	QPoint diff = e->pos() - mAnchor;
	 if (e->buttons() & Qt::LeftButton) {
		 mRotX += diff.y() * speedRot;
		 mRotY += diff.x() * speedRot;
		 print_cam();
		 updateGL();
	 }
	 if (e->buttons() & Qt::MiddleButton) { // Maya buttons
		 mCamPos.x += diff.x() * speedPan;
		 mCamPos.y -= diff.y() * speedPan;
         print_cam();
		 updateGL();
	 }
     if ( e->buttons() & Qt::RightButton ) {
         const float speed = 0.002f;
         mCamPos.z += speed * diff.y();
         print_cam();
         updateGL();
     }

	 mAnchor = e->pos();     
}

void GLWidget::mousePressEvent(QMouseEvent* e)
{
	mDownPos = mAnchor = e->pos();
}

void GLWidget::wheelEvent(QWheelEvent* e)
{
	const float speed = 0.002f;
	mCamPos.z += speed*e->delta();
    print_cam();
	updateGL();
}

void GLWidget::timerEvent(QTimerEvent* e)
{
	bool doRepaint = false;

	float speed = 0.005f;
	if (mMoveFast) speed *= 5.;
	
	if (mMoveState[MoveLeft])  { mCamPos.x += speed; doRepaint = true; }
	if (mMoveState[MoveRight]) { mCamPos.x -= speed; doRepaint = true; }
	if (mMoveState[MoveUp])    { mCamPos.y -= speed; doRepaint = true; }
	if (mMoveState[MoveDown])  { mCamPos.y += speed; doRepaint = true; }
	if (mMoveState[MoveOut])   { mCamPos.z -= speed; doRepaint = true; }
	if (mMoveState[MoveIn])    { mCamPos.z += speed; doRepaint = true; }
	if (doRepaint) 
		updateGL();
}

void GLWidget::print_cam() {
    if ( 0 ) {
        printf( "mCamPos = ( %g, %g, %g )", mCamPos.x, mCamPos.y, mCamPos.z );
        printf( ", mRot = ( %g, %g )\n", mRotX, mRotY );
    }
}

void GLWidget::setViewport(const Vec3i& gridsize) {
	if (mGridsize.x != gridsize.x ||
		mGridsize.y != gridsize.y ||
		mGridsize.z != gridsize.z) {
		if (mPlane < 0) {
			mPlane = gridsize[mPlaneDim] / 2;
		} else {
			Real fac = (Real)gridsize[mPlaneDim] / (Real)mGridsize[mPlaneDim];
			mPlane = (int)(fac * mPlane);
		}
		mGridsize = gridsize;
		Q_EMIT painterEvent(Painter::EventSetMax, mGridsize[mPlaneDim]);
		Q_EMIT painterEvent(Painter::EventSetPlane, mPlane);
	}
}

void GLWidget::keyPressEvent(QKeyEvent* e)
{
	if(!keyProcess(e->key(), e->modifiers(), true)) 
		QGLWidget::keyPressEvent(e);
	else 
		updateGL();
}

void GLWidget::keyReleaseEvent(QKeyEvent* e)
{
	if(!keyProcess(e->key(), e->modifiers(), false))
		QGLWidget::keyReleaseEvent(e);
	else
		updateGL();
}

bool GLWidget::keyProcess(int key, int modifier, bool down) 
{
	bool shift = (modifier & Qt::ShiftModifier);
	bool alt   = (modifier & Qt::AltModifier); 
	bool ctrl  = (modifier & Qt::ControlModifier); 
	if      (key == Qt::Key_A) { mMoveState[MoveLeft]  = down; mMoveFast = shift; }
	else if (key == Qt::Key_D) { mMoveState[MoveRight] = down; mMoveFast = shift; }
	else if (key == Qt::Key_W) { mMoveState[MoveUp]    = down; mMoveFast = shift; }
	else if (key == Qt::Key_S) { mMoveState[MoveDown]  = down; mMoveFast = shift; }
	else if (key == Qt::Key_Q) { mMoveState[MoveIn]    = down; mMoveFast = shift; }
	else if (key == Qt::Key_E) { mMoveState[MoveOut]   = down; mMoveFast = shift; }
	else if (down) 
	{
		// only press events
		// note Key_P and Key_L used for play/step in mainwindow.cpp
		if      (key == Qt::Key_Z)                  { /* next "solver" info sometime? */ }
		else if (key == Qt::Key_G)                  { Q_EMIT painterEvent(Painter::EventToggleGridDisplay); }
		// data grids, first int
		else if (key == Qt::Key_X && shift)         { /* int display mdoes, not yet used */ }
		else if (key == Qt::Key_X)                  { Q_EMIT painterEvent(Painter::EventNextInt);  updatePlane(mPlane); }

		else if ( key == Qt::Key_C && ctrl )
			exit(0);

		// real
		else if (key == Qt::Key_C && shift)         { Q_EMIT painterEvent(Painter::EventNextRealDisplayMode); /* real display modes */ }
		else if (key == Qt::Key_C)                  { Q_EMIT painterEvent(Painter::EventNextReal); updatePlane(mPlane); } 
		else if (shift && ((key == Qt::Key_Less) ||      
			    (key == Qt::Key_Comma) ) )          { Q_EMIT painterEvent(Painter::EventScaleRealDownSm); }
		else if (shift && ((key == Qt::Key_Greater) ||
			    (key == Qt::Key_Period) ) )         { Q_EMIT painterEvent(Painter::EventScaleRealUpSm); }
		else if ((key == Qt::Key_Less) ||      
			    (key == Qt::Key_Comma) )            { Q_EMIT painterEvent(Painter::EventScaleRealDown); }
		else if ((key == Qt::Key_Greater) ||
			    (key == Qt::Key_Period) )           { Q_EMIT painterEvent(Painter::EventScaleRealUp); }

		// vec3 grids, scaling can be used with two key combinations (the second one is for international keyboards)
		else if (key == Qt::Key_V && shift)         { Q_EMIT painterEvent(Painter::EventNextVecDisplayMode); }
		else if (key == Qt::Key_V)                  { Q_EMIT painterEvent(Painter::EventNextVec);  updatePlane(mPlane); }
		// grid scaling
		else if (key == Qt::Key_Colon)              { Q_EMIT painterEvent(Painter::EventScaleVecDownSm); }
		else if (key == Qt::Key_QuoteDbl)           { Q_EMIT painterEvent(Painter::EventScaleVecUpSm); }
		else if (key == Qt::Key_Semicolon)          { Q_EMIT painterEvent(Painter::EventScaleVecDown); }
		else if (key == Qt::Key_Apostrophe)         { Q_EMIT painterEvent(Painter::EventScaleVecUp); }

		// particles
		else if (key == Qt::Key_B && shift)         { Q_EMIT painterEvent(Painter::EventNextParticleDisplayMode); }
		else if (key == Qt::Key_B && alt)           { Q_EMIT painterEvent(Painter::EventNextSystem); }
		else if (key == Qt::Key_B)                  { Q_EMIT painterEvent(Painter::EventToggleParticles); }

		else if((key == Qt::Key_ParenLeft) ||       // a bit ugly, but for some reason parentheses dont work in some cases... fall back with dual assignment
			    (key == Qt::Key_9) )                { Q_EMIT painterEvent(Painter::EventScalePdataDown); }
		else if((key == Qt::Key_ParenRight) ||
			    (key == Qt::Key_0) )                { Q_EMIT painterEvent(Painter::EventScalePdataUp);   }

		// mesh display
		else if (key == Qt::Key_M && shift)           Q_EMIT painterEvent(Painter::EventMeshMode);
		else if (key == Qt::Key_BraceLeft )           { Q_EMIT painterEvent(Painter::EventScaleMeshDown); }
		else if (key == Qt::Key_BracketLeft)          { Q_EMIT painterEvent(Painter::EventScaleMeshDown); }
		else if (key == Qt::Key_BraceRight)           { Q_EMIT painterEvent(Painter::EventScaleMeshUp); }
		else if (key == Qt::Key_BracketRight)         { Q_EMIT painterEvent(Painter::EventScaleMeshUp); }
		// special mesh display modes
		else if (key == Qt::Key_M && alt)             Q_EMIT painterEvent(Painter::EventMeshColorMode);
		else if (key == Qt::Key_M && ctrl)            Q_EMIT painterEvent(Painter::EventToggleBackgroundMesh); 
		else if (key == Qt::Key_M)                    Q_EMIT painterEvent(Painter::EventNextMesh);
		
		// switch display planes
		else if ( (key == Qt::Key_Asterisk) || (key == Qt::Key_8) ) {
			mPlaneDim = (mPlaneDim+1) % 3;            
			Q_EMIT painterEvent(Painter::EventSetDim, mPlaneDim);
			Q_EMIT painterEvent(Painter::EventSetMax, mGridsize[mPlaneDim]);
		} 
		// move plane (+shift is faster)
		else if (shift && ((key == Qt::Key_Plus)  || (key == Qt::Key_Equal)) ) { 
			updatePlane(mPlane + 10);
		} else if (shift && ((key == Qt::Key_Minus) || (key == Qt::Key_Underscore)) ) { 
			updatePlane(mPlane - 10);
		} else if (key == Qt::Key_Plus || key == Qt::Key_Equal) { 
			updatePlane(mPlane + 1);
		} else if (key == Qt::Key_Minus) { 
			updatePlane(mPlane - 1);
		}
		else if ( key == Qt::Key_K) {
			QString filename = QString("scr_%1.png").arg(QString::number(mScreenshotNumber), 3, QChar('0'));
			screenshot(filename);
			mScreenshotNumber++;
		}
		
		else return false;
	}
	else return false;
	return true;
}

void GLWidget::screenshot(QString file) {
    paintGL();
	updateGL();
	glFlush();
	grabFrameBuffer().save(file, nullptr, -1);
}

void GLWidget::updatePlane(int plane) {
	mPlane = clamp(plane, 0, mGridsize[mPlaneDim]);
	Q_EMIT painterEvent(Painter::EventSetPlane, mPlane);
}



}
