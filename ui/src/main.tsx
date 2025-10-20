import { BrowserRouter } from 'react-router-dom'
import { createRoot } from 'react-dom/client'

import { Atlas } from './components/atlas';

import './style.css'
import 'bootstrap/dist/css/bootstrap.css';
import 'bootstrap/dist/js/bootstrap.bundle.js';
import 'bootstrap-icons/font/bootstrap-icons.css';



createRoot(document.getElementById('root')!).render(
  <BrowserRouter basename='/ui'>
    <Atlas />
  </BrowserRouter>,
)
