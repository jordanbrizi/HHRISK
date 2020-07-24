// Modules to control application life and create native browser window
const {app, BrowserWindow, Menu, ipcMain, dialog} = require('electron')
const { fstat } = require('fs')
const path = require('path')
const resultsPath = app.getAppPath() + '\\bin\\Results\\'
const appPath = app.getAppPath() + '\\'

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const createWindow = () => {
	const win = new BrowserWindow({
		width: 360,
		height: 640,
		backgroundColor: '#000',
		resizable: false,
		frame: false,
		show: true,
		icon: __dirname + '/favicon.ico',
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
		icon: __dirname + '/favicon.ico',
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
	
	// win.openDevTools()
	// winResults.openDevTools()

	ipcMain.on('resultados', () => {
		if (winResults.isVisible() == true) {
			winResults.hide()			
		} else {
			winResults.show()
		}
	})

	ipcMain.on('sair', () => app.quit())

	// -------------------------------------------------------------------------
	// -------------------------------------------------------------------------

	const obterJsons = () => {
		const fs = require('fs')
		const files = []
		fs.readdirSync(resultsPath).forEach(arquivo => {
			files.push(`${arquivo}`)
		})

		const jsons = files.filter(a => a.includes('.json'))

		return jsons
	}

	const Dados = () => {
		const hhr = ['Doses, HQ and CR.json',
			'Summary Aggregate Risk.json',
			'Extensive Aggregate Risk.json',
			'Summary Cumulative Risk.json',
			'Extensive Cumulative Risk.json',
			'Complementary Analyzes.json'
		] // Human Health Risk
		const er = ['Combined.json', 'Individual.json'] // Ecological Risk
		const rr = ['Radiological Risk.json'] // Radiological Risk

		const arquivos = {
			hhr: obterJsons().filter(a => hhr.includes(a)),
			er: obterJsons().filter(a => er.includes(a)),
			rr: obterJsons().filter(a => rr.includes(a))
		}

		return {
			quantidade: () => arquivos.length,
			arquivos: () => arquivos,
			pegar: arquivo => require(path.resolve(resultsPath + arquivo))
		}
	}

	ipcMain.on('gerarOds', (event, arg) => {
		const xlsx = require('xlsx')
		const planilhas = []
		jsons = Dados().arquivos()[arg.grupo] // PEGA OS JSONS DO GRUPO t
		jsons.forEach(json => {
			arquivo = Dados().pegar(json)
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
		obterJsons().forEach(json => fs.unlinkSync(resultsPath + json))
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
	})
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

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