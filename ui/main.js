// Modules to control application life and create native browser window
const {app, BrowserWindow, Menu, ipcMain, dialog} = require('electron')
const { fstat } = require('fs')
const createWindow = () => {
	const win = new BrowserWindow({
		width: 360,
		height: 640,
		backgroundColor: '#000',
		resizable: false,
		frame: false,
		show: false,
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
		show: true,
		webPreferences: {
			nodeIntegration: true
		}
	})

	winResults.loadURL(`file://${__dirname}/results.html`)
	win.loadURL(`file://${__dirname}/index.html`)

	// win.once('ready-to-show', () => {
	// 	win.show()
	// })

	Menu.setApplicationMenu(null)
	
	//win.openDevTools()
	winResults.openDevTools()

	ipcMain.on('resultados', () => {
		if (winResults.isVisible() == true) {
			winResults.hide()			
		} else {
			winResults.show()
		}
	})

	ipcMain.on('sair', () => app.quit())

	ipcMain.on('salvar_planilha', (event, arg) => {
		let path = require('path')
		let options = {
			title: "Selecionar Pasta",
			defaultPath: path.resolve(require('os').homedir()),
			properties: ['openDirectory']
		}
		dialog.showOpenDialog(options).then((response) => {
			if(response.canceled === false) {
				const fs = require('fs')
				arg.sheets.forEach(arquivo => {
					file_name = path.resolve(
						response.filePaths + arquivo.replace(arg.path, '/')
					)
					fs.rename(arquivo, file_name, err => {
						if(err) throw err
					})
				})
				require('child_process')
					.exec(`start "" "${response.filePaths}"`)
			}
		}).catch(err => {
			console.log(err)
		})
	})
}

app.on('ready', createWindow)
app.on('window-all-closed', () => {
	if (process.platform !== 'darwin') app.quit()
})
app.on('activate', function () {
	if (BrowserWindow.getAllWindows().length === 0) {
		createWindow()
	}
})