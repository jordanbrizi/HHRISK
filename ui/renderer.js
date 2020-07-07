const { ipcRenderer } = require('electron')
const pkg = () => require('./package')
const tabela = () => require('./tabela')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale

const lang = (attr = defaultLang) => {
	if(attr != null) {
		return require(`./lang/${attr}`)
	}
	else {
		return require(`./lang/${defaultLang}`)	
	}
}
module.exports = {pkg, lang, tabela, ipcRenderer}