// Modules to control application life and create native browser window
const {app, BrowserWindow, Menu, ipcMain} = require('electron')
const createWindow = () => {
	const win = new BrowserWindow({
		width: 360,
		height: 640,
		backgroundColor: '#000',
		resizable: false,
		frame: false,
		webPreferences: {
			nodeIntegration: true
		}
	})
	const winResults = new BrowserWindow({
		width: 1200,
		height: 675,
		backgroundColor: '#000',
		resizable: true,
		frame: false,
		titleBarStyle: 'hidden',
		show: false,
		webPreferences: {
			nodeIntegration: true
		}
	})

	winResults.loadURL(`file://${__dirname}/results.html`)
	win.loadURL(`file://${__dirname}/index.html`)

	win.once('ready-to-show', () => {
		win.show()
	})

	Menu.setApplicationMenu(null)
	win.openDevTools()
	winResults.openDevTools()

	ipcMain.on('resultados', () => {
		if (winResults.isVisible() == true) {
			winResults.hide()			
		} else {
			winResults.show()
		}
	})
}

app.on('ready', createWindow)
app.on('window-all-closed', function () {
	if (process.platform !== 'darwin') {
		app.quit()
	}
})
app.on('activate', function () {
	if (BrowserWindow.getAllWindows().length === 0) {
		createWindow()
	}
})