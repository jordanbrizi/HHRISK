<<<<<<< HEAD
const {app, BrowserWindow, Menu, ipcMain, dialog} = require('electron')
const { fstat } = require('fs')
const path = require('path')
const resultsPath = app.getAppPath() + '\\bin\\Results\\'
const appPath = app.getAppPath() + '\\'

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

=======
// Modules to control application life and create native browser window
const {app, BrowserWindow, Menu, ipcMain} = require('electron')
>>>>>>> 39d9aa5686bc26510bb8baebdb57f262b29ce598
const createWindow = () => {
	const win = new BrowserWindow({
		width: 360,
		height: 640,
		backgroundColor: '#000',
		resizable: false,
		frame: false,
<<<<<<< HEAD
		show: true,
		icon: __dirname + '/favicon.ico',
=======
>>>>>>> 39d9aa5686bc26510bb8baebdb57f262b29ce598
		webPreferences: {
			nodeIntegration: true
		}
	})
	const winWarnings = new BrowserWindow({
		width: 360,
		height: 640,
		backgroundColor: '#FFF',
		resizable: false,
		show: false,
		icon: __dirname + '/favicon.ico',
		webPreferences: {
			nodeIntegration: true
		}
	})
	const winInformation = new BrowserWindow({
		width: 360,
		height: 640,
		backgroundColor: '#FFF',
		resizable: false,
		show: false,
		icon: __dirname + '/favicon.ico',
		webPreferences: {
			nodeIntegration: true
		}
	})

	win.loadURL(`file://${__dirname}/index.html`)
	winWarnings.loadURL(resultsPath + 'WARNINGS.txt')
	winInformation.loadURL(resultsPath + 'Information.txt')

	win.once('ready-to-show', () => {
		win.show()
	})
	winWarnings.on('closed', () => {
		winWarnings = null
	})

	Menu.setApplicationMenu(null)
	
<<<<<<< HEAD
	// win.openDevTools()

	ipcMain.on('sair', () => app.quit())

	// -------------------------------------------------------------------------
	// -------------------------------------------------------------------------
		
	const Obter = () => {
		const files = []
		const fs = require('fs')
		fs.readdirSync(resultsPath).forEach(arquivo => {
			files.push(`${arquivo}`)
		})

		return {
			jsons: files.filter(a => a.includes('.json')),
			txts: files.filter(a => a.includes('.txt')),
		}
	}

	const Resultados = {
		jsons: Obter().jsons,
		txts: Obter().txts,
		quantidade: tipo => Resultados[tipo].length,
		arquivos: () => arquivos,
		pegar: arquivo => require(path.resolve(resultsPath + arquivo))
	}

	ipcMain.on('gerarOds', (event, arg) => {
		const xlsx = require('xlsx')
		const planilhas = []
		jsons = Resultados.jsons
		jsons.forEach(json => {
			arquivo = Resultados.pegar(json)
			chaves = Object.keys(arquivo)
			const wb = xlsx.utils.book_new()
			chaves.forEach(chave => {
				chaveNew = chave.substring(0, 28) + '...' //C/ ATÉ 31 CARACTERES
				keys = Object.keys(arquivo[chave][0])
				header = [{chave: chave}]
				ws = xlsx.utils.json_to_sheet(header, { skipHeader: true })
				xlsx.utils.sheet_add_json(ws, arquivo[chave], { origin: "A2" })
				const merge = [{ s: { r: 0, c: 0 }, e: { r: 0, c: (keys.length -1) } }]
				ws["!merges"] = merge
				xlsx.utils.book_append_sheet(
					wb,
					ws,
					chaveNew
				)
			})
			let sheetPath = app.getPath('temp')
			let sheetName = `\\${json.replace('.json', '')}.ods`
			xlsx.writeFile(wb, sheetPath+sheetName)
			planilhas.push(sheetName)
		})
		let options = {
			title: "Selecionar Pasta",
			defaultPath: app.getPath('documents'),
			properties: ['openDirectory']
		}

	// ABRIR O DIÁLOGO DE SELEÇÃO DE PASTA
		dialog.showOpenDialog(options).then((response) => {
			if (response.canceled === false) {
				const fs = require('fs')
				planilhas.forEach(sheet => {
					oldPath = path.resolve(app.getPath('temp') + sheet)
					newPath = path.resolve(response.filePaths + sheet)
					fs.rename(oldPath, newPath, err => {
						if (err) throw err
					})
				})
				require('child_process')
					.exec(`start "" "${response.filePaths}"`)
			}
		}).catch(err => {
			console.log(err)
		})
	})

	ipcMain.on('clearResults', (event, arg) => {
		const fs = require('fs')
		Resultados.jsons.forEach(json => fs.unlinkSync(resultsPath + json))
		return
	})

	ipcMain.on('execute', (event, arg) => {
		var child = require('child_process')
		var path = require('path')
		var hhrisk_exe = appPath + 'bin\\HERisk.exe'
		child.exec(`cd "${appPath}bin" & cmd /K ${hhrisk_exe}`, (err, data, t) => {
			if (err) {
				console.error(err)
				return
			}
		})
		event.sender.send('executed', true)
=======
	win.openDevTools()
	winResults.openDevTools()

	ipcMain.on('resultados', () => {
		if (winResults.isVisible() == true) {
			winResults.hide()			
		} else {
			winResults.show()
		}
	})
	ipcMain.on('sair', () => {
		app.quit()
>>>>>>> 39d9aa5686bc26510bb8baebdb57f262b29ce598
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
app.setAboutPanelOptions({
	applicationName: "HERisk",
	applicationVersion: app.getVersion(),
	copyright: "Todos os direitos reservados",
	version: app.getVersion(),
	iconPath: appPath + 'ui\\favicon.ico'
})