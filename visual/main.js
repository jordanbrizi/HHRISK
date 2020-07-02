// Modules to control application life and create native browser window
const {app, BrowserWindow, Menu} = require('electron')
console.log('😍')
function createWindow () {
	// Create the browser window.
	const win = new BrowserWindow({
		width: 360,
		height: 640,
		title: 'HHRISK',
		show: false,
		backgroundColor: '#283E51',
		resizable: false,
		titleBarStyle: 'hidden',
		frame: false,
		webPreferences: {
			nodeIntegration: true
		}
	})

	win.once('ready-to-show', () => {
		win.show()
	})
	// and load the index.html of the app.
	win.loadFile('index.html')
	Menu.setApplicationMenu(null)
	// Open the DevTools.
	win.webContents.openDevTools()
}

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', createWindow)

// Quit when all windows are closed.
app.on('window-all-closed', function () {
	// On OS X it is common for applications and their menu bar
	// to stay active until the user quits explicitly with Cmd + Q
	if (process.platform !== 'darwin') {
		app.quit()
	}
})

app.on('activate', function () {
	// On OS X it's common to re-create a window in the app when the
	// dock icon is clicked and there are no other windows open.
	if (BrowserWindow.getAllWindows().length === 0) {
		createWindow()
	}
})

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.