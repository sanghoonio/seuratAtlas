import { EmbeddingAtlas } from 'embedding-atlas/react';
import * as vg from '@uwdata/vgplot'

export function Atlas() {

  const coordinator = new vg.Coordinator(vg.socketConnector('ws://localhost:3000/' as any));
  
  return (
    <div className='vh-100'>
      <EmbeddingAtlas
        coordinator={coordinator}
        data={{
          table: 'seurat_obj',
          id: 'cell_id',
          projection: { x: 'projection_x', y: 'projection_y' },
        }}
      />
    </div>
    
  )
}