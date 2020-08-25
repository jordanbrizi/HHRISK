const { ipcRenderer, remote, shell } = require('electron')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale
const mainPath = arg => remote.app.getAppPath() + '\\' + arg

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

module.exports = { defaultLang, ipcRenderer, shell, remote }