/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * QT main window
 *
 ******************************************************************************/

#include "mainwindow.h"
#include "qtmain.h"

#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QtOpenGL>
#include <sstream>
#include "meshpainter.h"
#include "particlepainter.h"

using namespace std;

namespace Manta {

MainWnd::MainWnd() : QMainWindow( 0 ), mPaused( true ), mRequestPause( false ), mRequestClose( false ), mStep( 0 ),
                     mKbwScene( 0 ), mKbwView( 0 ), mKbwPixmap( 0 ), mMenuBar( 0 ) {
    // Frame info label
    mInfo = new QLabel;

    QFont font( "Consolas", 12 );
    QApplication::setFont( font, "QLabel" );
    // mInfo->setFont( font ); // doesn't change size

    setStep( 0, 0. );

    // register GL widget
    mGlWidget = new GLWidget();
    setCentralWidget( mGlWidget );
    connect( mGlWidget, SIGNAL( clickLine( QPoint, float, float, float, float, float, float ) ), SLOT( clickLine( QPoint, float, float, float, float, float, float ) ) );

    // register grid painters
    mPainterLayout = new QVBoxLayout;
    mPainterLayout->setAlignment( Qt::AlignTop );
    mPainterLayout->addWidget( mInfo );
    GridPainter<int> *intPainter = new GridPainter<int>( NULL, this );
    mPainter.push_back( new GridPainter<Real>( (FlagGrid **)intPainter->getGridPtr(), this ) );
    if ( 1 ) { // 1(default); different painting order
        mPainter.push_back( new GridPainter<Vec3>( NULL, this ) );
        mPainter.push_back( intPainter );
    } else {
        mPainter.push_back( intPainter );
        mPainter.push_back( new GridPainter<Vec3>( NULL, this ) );
    }	
	mPainter.push_back(new ParticlePainter(intPainter, this));
	MeshPainter* ptr = new MeshPainter(this);
	mPainter.push_back(ptr);    
	connect(this, SIGNAL(setBackgroundMesh(Mesh*)), ptr, SLOT(setBackgroundMesh(Mesh*)));

	for (int i=0; i<(int)mPainter.size(); i++) {
		connect(mGlWidget, SIGNAL(paintSub()), mPainter[i], SLOT(paint()));
		connect(mGlWidget, SIGNAL(painterEvent(int, int)), mPainter[i], SLOT(doEvent(int, int)));
		connect(this, SIGNAL(painterEvent(int, int)), mPainter[i], SLOT(doEvent(int, int)));
		connect(mPainter[i], SIGNAL(setViewport(const Vec3i&)), mGlWidget, SLOT(setViewport(const Vec3i&)));
		mPainter[i]->attachWidget(mPainterLayout);
	}
	
	// docking widget for painters
	QDockWidget* painterDock = new QDockWidget("Info", this);
	QWidget* painterProxy = new QWidget;
	painterProxy->setLayout(mPainterLayout);    
	painterDock->setWidget(painterProxy);
	painterDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	addDockWidget(Qt::RightDockWidgetArea, painterDock);
	 
	// Top toolbar
	QToolBar* toolbar = addToolBar("Control");
	toolbar->setAllowedAreas(Qt::TopToolBarArea);
	toolbar->setMovable(false);    
	mAcPlay = toolbar->addAction(QIcon(":/play.png"),"Play");
	mAcPlay->setStatusTip("Continue simulation");
	connect(mAcPlay, SIGNAL(triggered()), SLOT(play()));
	mAcPause = toolbar->addAction(QIcon(":/pause.png"),"Pause");
	mAcPause->setStatusTip("Pause simulation");
	connect(mAcPause, SIGNAL(triggered()), SLOT(pause()));
	Q_EMIT play();
	
	// build menu
	/*QAction* a = new QAction(this);
	a->setText( "Quit" );
	connect(a, SIGNAL(triggered()), SLOT(close()) );
	mMenuBar = menuBar()->addMenu( "File" );
	mMenuBar->addAction( a ); */

	// keyboard info window, show on demand
    mKbwScene = new QGraphicsScene(); 
    mKbwView = new QGraphicsView(mKbwScene);
    mKbwPixmap = new QGraphicsPixmapItem(QPixmap(":/keyboard.png"));
    mKbwScene->addItem(mKbwPixmap);
	mKbwView->hide(); 

	mAcHelp = toolbar->addAction(QIcon(":/help.png"),"Help");
	mAcHelp->setStatusTip("Help");
	connect(mAcHelp, SIGNAL(triggered()), SLOT(showHelp()));
	
	// start...
	mGlWidget->setFocus();
	this->raise();
	this->activateWindow();

	// move gui window to upper left corner and resize window to screen size, enable on demand...
	if(false) {
		QRect rc = frameGeometry();
		QRect rcDesktop = QApplication::desktop()->frameGeometry();
		rc.setLeft(rcDesktop.left());
		rc.setTop(rcDesktop.top());
		rc.setRight(rcDesktop.right());
		rc.setBottom(rcDesktop.bottom());
		move(rc.topLeft());
		resize(rc.size());
	}

	// uncomment to start  paused
	//Q_EMIT pause();
}

void MainWnd::clickLine(QPoint pos, float p0, float p1,float p2, float q0, float q1, float q2) {
	string msg;
	for (int i=mPainter.size()-1; i>=0; i--) {
		msg += mPainter[i]->clickLine(Vec3(p0,p1,p2),Vec3(q0,q1,q2));
	}
	if (!msg.empty())
		QToolTip::showText(pos, QString(msg.c_str()));
}

void MainWnd::addControl(void* ctrl) {
	CustomControl* control = (CustomControl*) ctrl;
	mCtrls.push_back(control);
	control->init(mPainterLayout);
}

void MainWnd::setStep(int f, float time) {
	// weird crash on QAccessible::updateAccessibility() if f or time vary
	//f = 1; time = 1;

    if ( 0 ) { 
        std::stringstream s;
        s << "Simulation frame " << f << "\nTime " << time;
        auto s1 = s.str();
        auto s2 = s1.c_str();
        mInfo->setText( s2 );
    } else {
		static char s[100];
		//sprintf( s, "Simulation frame %d\nTime %g", f, time );
        //sprintf( s, "Time %g", time );
        sprintf( s, "Time %d(count=%0.1f)", f, time );
		mInfo->setText( s );
    }
}

void MainWnd::setPauseStatus(bool v)
{
	mPaused = v; 
}

bool MainWnd::event(QEvent* e) {
	if (e->type() == (QEvent::Type)EventGuiShow) {
		if (!mRequestClose) {
			//this->show();

			this->showMaximized();
            //this->setWindowState( Qt::WindowMaximized );

			//int screen_height = QApplication::desktop()->screenGeometry().height();
            //int screen_width = QApplication::desktop()->screenGeometry().width();
            //this->setMinimumSize( screen_width - 10, screen_height - 10 );

            Q_EMIT painterEvent( Painter::UpdateFull );
            mGlWidget->updateGL();

			// if happens before the paint, then I get "QOpenGLContext::swapBuffers() called with non-exposed window, behavior is undefined"
			this->showMinimized();
		}
		Q_EMIT wakeMain();
		return true;
	}
	else if (e->type() == (QEvent::Type)EventFullUpdate) {
		if (!mRequestClose) {
			Q_EMIT painterEvent(Painter::UpdateFull);
			mGlWidget->updateGL();
		}
		Q_EMIT wakeMain();
		return true;
	}
	else if (e->type() == (QEvent::Type)EventStepUpdate) {        
		if (!mRequestClose) {
			if (mRequestPause) {
				Q_EMIT painterEvent(Painter::UpdateFull);
				mGlWidget->updateGL();
			} else {
				Q_EMIT painterEvent(Painter::UpdateStep);
				// redraw not necessary? old: mGlWidget->updateGL();
			}
		}
		Q_EMIT wakeMain();
		return true;
	}
	else if (e->type() == (QEvent::Type)EventFinalUpdate) {        
		if (!mRequestClose) {
			Q_EMIT painterEvent(Painter::UpdateFull);
			mGlWidget->updateGL();
		}
		mRequestClose = true;
		Q_EMIT wakeMain();
		return true;
	}
	else if (e->type() == (QEvent::Type)EventSet2DCam) {        
		mGlWidget->setCamPos( Vec3(0, 0, -1.3) );
		return true;
	}
	else if (e->type() == (QEvent::Type)EventInstantKill) {        
		Q_EMIT killMain();
		Q_EMIT exitApp();
		return true;
	}

	// update button states for pause events
	if( (mRequestPause) && (!mAcPlay->isEnabled()) ) {
			mAcPlay->setEnabled(true);
			mAcPause->setEnabled(false);    
	}
	if( (mRequestPause) && (!mAcPlay->isEnabled()) ) {
			mAcPlay->setEnabled(true);
			mAcPause->setEnabled(false);    
	}
	
	return QMainWindow::event(e);
}

// there's also GLWidget::keyProcess()
void MainWnd::keyPressEvent(QKeyEvent* e) {
	if ( e->key() == Qt::Key_Escape ) {
		exit(0); // zl
		mRequestClose = true;
		Q_EMIT killMain();
		this->close();
	} else if (0&& e->key() == Qt::Key_Space) {
		if (mRequestClose) {
			Q_EMIT killMain();
			this->close();
		} else {
			Q_EMIT painterEvent(mPaused ? Painter::UpdateFull : Painter::UpdateRequest);
			mGlWidget->updateGL();
		}
	} else if (e->key() == Qt::Key_P || e->key() == Qt::Key_Space) {
		if (mRequestClose) {
			Q_EMIT killMain();
			this->close();
		} else if (mRequestPause)
			Q_EMIT play();
		else    
			Q_EMIT pause();
	} else if (e->key() == Qt::Key_L || e->key() == Qt::Key_Right) {
		if (mRequestClose) {
			Q_EMIT killMain();
			this->close();
		} else if (mRequestPause) {
			mRequestPause = false;
			mStep = (e->modifiers() & Qt::ShiftModifier) ? 1 : 2;                
		} else
			Q_EMIT pause();
	} else if (e->key() == Qt::Key_H) {
		Q_EMIT showHelp();
	} else { 
		mGlWidget->keyPressEvent(e); // let gl widget take care of keyboard shortcuts
		//QMainWindow::keyPressEvent(e);
	}
}
void MainWnd::keyReleaseEvent(QKeyEvent* e)
{
	mGlWidget->keyReleaseEvent(e);
}

void MainWnd::pause() {
	mRequestPause = true;
	// dont call: mAcPlay/mAcPause ->setEnabled(true) here; wrong thread if called from python
}

void MainWnd::play() {
	mRequestPause = false;
	mAcPlay->setEnabled(false);
	mAcPause->setEnabled(true);    
}

void MainWnd::step() {
	mStep = 2;
	mRequestPause = false;
}

void MainWnd::showHelp() {
	mKbwView->show();
}

void MainWnd::setRealGridDisplay( int dm ) {
	Q_EMIT painterEvent(Painter::EventSetRealGridDisplayMode, dm); 
}

void MainWnd::setVec3GridDisplay( int dm ) {
	Q_EMIT painterEvent(Painter::EventSetVecGridDisplayMode, dm); 
}

void MainWnd::nextRealGrid() {
	Q_EMIT painterEvent(Painter::EventNextReal); 
}
void MainWnd::nextVec3Grid() {
	Q_EMIT painterEvent(Painter::EventNextVec); 
}
void MainWnd::nextMesh() {
	Q_EMIT painterEvent(Painter::EventNextMesh); 
}
void MainWnd::nextParts() {
	Q_EMIT painterEvent(Painter::EventNextSystem); 
}
void MainWnd::nextPdata() {
	Q_EMIT painterEvent(Painter::EventToggleParticles); 
}
void MainWnd::nextVec3Display() {
	Q_EMIT painterEvent(Painter::EventNextVecDisplayMode); 
}
void MainWnd::nextPartDisplay() {
	Q_EMIT painterEvent(Painter::EventNextParticleDisplayMode); 
}
void MainWnd::nextMeshDisplay() {
	Q_EMIT painterEvent(Painter::EventMeshMode); 
}
void MainWnd::toggleHideGrids() {
	Q_EMIT painterEvent(Painter::EventToggleGridDisplay); 
}
void MainWnd::setCamPos(float x, float y, float z) {
	mGlWidget->setCamPos( Vec3(x, y, z) );
}
void MainWnd::setCamRot(float x, float y, float z) {
	mGlWidget->setCamRot( Vec3(x, y, z) );
}

void MainWnd::setPlane( int plane ) {
    mGlWidget->setPlane( plane );
}

void MainWnd::windowSize(int w, int h) {
	mGlWidget->setMinimumSize( w,h );
	mGlWidget->setMaximumSize( w,h );
	mGlWidget->resize( w,h );
}

void MainWnd::showEvent_( QShowEvent *event ) {
    QMainWindow::showEvent( event );
	static int count = 0;
	// zl first time it's minimized on EventGuiShow
	// no good: opengl sets the point size based on the window size which needs to be maximized; disabling the function for now
    if ( 0&& ++count == 2 )
        this->showMaximized();
}

MainWnd::~MainWnd() {
}

void MainWnd::screenshot(QString file) {
	mGlWidget->screenshot(file);
}


}
