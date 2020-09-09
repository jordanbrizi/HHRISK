const { app, BrowserWindow, Menu, ipcMain, dialog, Notification } = require('electron')
const path = require('path')
const resultsPath = app.getAppPath() + '\\bin\\Results\\'
const appPath = app.getAppPath() + '\\'

// OBTER OS ARQUIVOS EM JSON E TXT DA PASTA RESULTS	
const Obter = () => require('fs').readdirSync(resultsPath)
// FILTRAR OS RESULTADOS EM JSON OU TXT
const Resultados = tipo => {
	switch (tipo) {
		case 'txts':
			return Obter().filter(a => a.includes('.txt'))
			break

		default:
			return Obter().filter(a => a.includes('.json'))
			break
	}
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const createWindow = () => {
	app.allowRendererProcessReuse = false
	setInterval(() => { Resultados() }, 1000)
	const win = new BrowserWindow({
		width: 360,
		height: 640,
		backgroundColor: '#23272A',
		resizable: false,
		frame: false,
		show: true,
		icon: __dirname + '/favicon.ico',
		webPreferences: {
			nodeIntegration: true
		}
	})	

	win.loadURL(`file://${__dirname}/index.html`)
	
	win.once('ready-to-show', () => win.show())
	
	Menu.setApplicationMenu(null)
	
	// win.openDevTools()
	
	ipcMain.on('how_to', () => {
		const winHowto = new BrowserWindow({
			width: 1006,
			height: 595,
			backgroundColor: '#23272A',
			resizable: false,
			show: false,
			icon: __dirname + '/favicon.ico',
			webPreferences: {
				nodeIntegration: true
			}
		})
		winHowto.loadURL('https://www.youtube.com/embed/w3-UVmYQrmM')
		winHowto.show()
	})
	ipcMain.on('guide', () => {
		const winGuide = new BrowserWindow({
			width: 1024,
			height: 640,
			backgroundColor: '#23272A',
			resizable: true,
			show: false,
			icon: __dirname + '/favicon.ico'
		})
		winGuide.loadURL(appPath + 'bin/HERisk.pdf')
		winGuide.show()
	})
	ipcMain.on('guideSave', () => {
		const options = {
			defaultPath: app.getPath('documents'),
			properties: ['openDirectory']
		}
		dialog.showOpenDialog(options).then((response) => {
			if (response.canceled === false) {
				const fs = require('fs')
				const old = path.resolve(appPath + 'bin\\HERisk.docx')
				const nFile = path.resolve(response.filePaths + '\\HERisk.docx')
				fs.copyFile(old, nFile, () => {
					require('child_process')
						.exec(`start "" "${response.filePaths}"`)
				})
			}
		}).catch(err => {
			event.sender.send('responseError', err)
		})
	})

	ipcMain.on('sair', () => app.quit())

	// -------------------------------------------------------------------------
	// -------------------------------------------------------------------------

	ipcMain.on('gerarOds', (event, arg) => {
		const messages = {
			default: [
				'Generating ODS files.',
				'Select Folder',
				'ODS files has been generated.',
				'Canceled by user.'
			],
			br: [
				'Gerando arquivos ODS.',
				'Selecionar pasta',
				'As planilhas ODS foram geradas.',
				'Cancelado pelo usuário.'
			]
		}
		event.sender.send('responseSuccess', messages[arg][0])
		const xlsx = require('xlsx')
		const planilhas = []
		const options = {
			title: messages[arg][1],
			defaultPath: app.getPath('documents'),
			properties: ['openDirectory']
		}
		jsons = Resultados()
		// console.log(jsons)
		jsons.forEach(json => {
			const arquivo = require(path.resolve(resultsPath + json))
			const chaves = Object.keys(arquivo)
			const wb = xlsx.utils.book_new()
			chaves.forEach(chave => {
				const chaveNew = chave.substring(0, 28) + '...' //C/ ATÉ 31 CARACTERES
				const keys = Object.keys(arquivo[chave][0])
				const header = [{ chave: chave }]
				const ws = xlsx.utils.json_to_sheet(header, { skipHeader: true })
				xlsx.utils.sheet_add_json(ws, arquivo[chave], { origin: "A2" })
				const merge =
					[{ s: { r: 0, c: 0 }, e: { r: 0, c: (keys.length - 1) } }]
				ws["!merges"] = merge
				xlsx.utils.book_append_sheet(wb, ws, chaveNew)
			})
			const sheetName = `\\${json.replace('.json', '')}.ods`
			xlsx.writeFile(wb, app.getPath('temp') + sheetName)
			planilhas.push(sheetName)
		})

		event.sender.send('responseSuccess', messages[arg][2])

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
			} else {
				event.sender.send('responseError', messages[arg][3])
			}
		}).catch(err => {
			console.log(err)
		})
	})

	ipcMain.on('execute', (event, arg) => {
		const messages = {
			default: [
				'HERisk is running...',
				'Cleaning the results folder...',
				'An error occurred in ### sheet. Please check the data entered and try again.',
				'A value of 0 was inserted in some field. Please check the entered values.',
				'Successfully executed.'
			],
			br: [
				'HERisk está sendo executado...',
				'Limpando os resultados anteriores...',
				'Ocorreu um erro na planilha ###. Por favor, verifique os dados inseridos e tente novamente.',
				'Um valor 0 ou nulo foi inserido em algum campo. Por favor, verifique os dados inseridos e tente novamente.',
				'Executado com sucesso.'
			]
		}
		event.sender.send('responseSuccess', messages[arg][0])
		const fs = require('fs')
		const child = require('child_process')
		const herisk_exe = appPath + 'bin\\HERisk.exe'

		// LIMPA A PASTA RESULTS
		event.sender.send('responseSuccess', messages[arg][1])
		if (Resultados().length > 0) {
			Resultados().forEach(json => fs.unlinkSync(resultsPath + json))
		}
		if (Resultados('txts').length > 0) {
			Resultados('txts').forEach(txt => fs.unlinkSync(resultsPath + txt))
		}
		child.exec(herisk_exe, { "cwd": appPath + "bin" }, (err, data, stderr) => {
			if (err) {
				event.sender.send('responseError', stderr)
				const prns = [
					"Concentration",
					"Datachemical",
					"Dataecological",
					"Dataexp",
					"Scenary"
				]
				const prn = prns.filter(a => stderr.includes(a))
				const erros = ["divide by zero"]
				const erro = erros.filter(a => stderr.includes(a))

				if (prn.length > 0) {
					event.sender.send('responseError',
						messages[arg][2].replace('###', prn))
				}
				if (erro.length > 0) {
					event.sender.send('responseError', messages[arg][3])
				}
			} else {
				event.sender.send('responseSuccess', messages[arg][4])
			}
		})
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