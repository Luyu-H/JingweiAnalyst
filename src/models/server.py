from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import Optional, Dict
import json
import asyncio
import uuid

from src.models.agent_with_tools_v3 import AgentWithTools
from src.tools.registry import TOOLS, auto_import_tools

auto_import_tools()

app = FastAPI(title="Scientific Agent API")

# CORS middleware for frontend access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify your frontend URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve static files (for visualizations)
app.mount("/static", StaticFiles(directory="./results"), name="static")

# Data configuration
DATA_LIST = [
    {
        "name": "nitrate_data_v1.1",
        "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_no3_sjtu_jingwei_v1_1.zarr",
        "description": "Nitrate concentration data from ocean reconstruction, version 1.1"
    },
    {
        "name": "oxygen_data_v1.1",
        "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr",
        "description": "Dissolved oxygen concentration data from ocean reconstruction, version 1.1"
    },
    {
        "name": "phosphate_data_v1.1", 
        "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_po4_sjtu_jingwei_v1_1.zarr",
        "description": "Phosphate concentration data from ocean reconstruction, version 1.1"
    },
    {
        "name": "silicate_data_v1.1", 
        "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_si_sjtu_jingwei_v1_1.zarr",
        "description": "Silicate concentration data from ocean reconstruction, version 1.1"
    },
    {
        "name": "temperature_data_v1.1", 
        "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_temp_sjtu_jingwei_v1_1.zarr",
        "description": "Sea temperature data from ocean reconstruction, version 1.1"
    },
    {
        "name": "salinity_data_v1.1", 
        "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_so_sjtu_jingwei_v1_1.zarr",
        "description": "Sea salinity data from ocean reconstruction, version 1.1"
    }
]

# Initialize agent (singleton)
API_KEY = "sk-edc62146735146e4a9904740057ace47"
agent = AgentWithTools(
    TOOLS, 
    DATA_LIST, 
    API_KEY,
    max_turns_per_session=10,
    max_code_iterations=3
)

# Active WebSocket connections
active_connections: Dict[str, WebSocket] = {}


class QueryRequest(BaseModel):
    """HTTP request model for queries."""
    query: str
    session_id: Optional[str] = None


class SessionRequest(BaseModel):
    """Request model for session management."""
    session_id: Optional[str] = None


@app.get("/")
async def root():
    """Root endpoint."""
    return {"message": "Scientific Agent API", "version": "1.0"}


@app.post("/api/session/create")
async def create_session(request: SessionRequest):
    """Create a new conversation session."""
    session_id = agent.create_session(request.session_id)
    return {"session_id": session_id, "status": "created"}


@app.post("/api/session/end")
async def end_session(request: SessionRequest):
    """End an existing session."""
    if not request.session_id:
        raise HTTPException(status_code=400, detail="session_id required")
    
    agent.end_session(request.session_id)
    return {"session_id": request.session_id, "status": "ended"}


@app.get("/api/session/{session_id}/history")
async def get_session_history(session_id: str):
    """Get conversation history for a session."""
    if session_id not in agent.sessions:
        raise HTTPException(status_code=404, detail="Session not found")
    
    session = agent.sessions[session_id]
    return {
        "session_id": session_id,
        "turn_count": session.turn_count,
        "max_turns": session.max_turns,
        "history": session.history
    }


@app.websocket("/ws/{client_id}")
async def websocket_endpoint(websocket: WebSocket, client_id: str):
    """
    WebSocket endpoint for real-time communication.
    
    Message format from client:
    {
        "type": "query",
        "query": "user question",
        "session_id": "optional session id"
    }
    
    Message format to client:
    {
        "type": "status",
        "stage": "planning|selection|execution|visualization|complete",
        "data": {...}
    }
    """
    await websocket.accept()
    active_connections[client_id] = websocket
    
    try:
        # Send welcome message
        await websocket.send_json({
            "type": "connected",
            "message": "Connected to Scientific Agent",
            "client_id": client_id
        })
        
        while True:
            # Receive message from client
            data = await websocket.receive_json()
            
            if data.get("type") == "query":
                query = data.get("query")
                session_id = data.get("session_id")
                
                if not query:
                    await websocket.send_json({
                        "type": "error",
                        "error": "Query cannot be empty"
                    })
                    continue
                
                # Status callback function
                async def status_callback(status_data: dict):
                    await websocket.send_json({
                        "type": "status",
                        **status_data
                    })
                
                try:
                    # Run agent with streaming updates
                    result = await agent.run_streaming(
                        query,
                        session_id=session_id,
                        status_callback=status_callback
                    )
                    
                    # Send completion message with full result
                    await websocket.send_json({
                        "type": "result",
                        "session_id": result.get("session_id"),
                        "success": True
                    })
                    
                except Exception as e:
                    await websocket.send_json({
                        "type": "error",
                        "error": str(e)
                    })
            
            elif data.get("type") == "ping":
                await websocket.send_json({"type": "pong"})
            
            else:
                await websocket.send_json({
                    "type": "error",
                    "error": f"Unknown message type: {data.get('type')}"
                })
    
    except WebSocketDisconnect:
        if client_id in active_connections:
            del active_connections[client_id]
        print(f"Client {client_id} disconnected")
    
    except Exception as e:
        print(f"WebSocket error for client {client_id}: {e}")
        if client_id in active_connections:
            del active_connections[client_id]


@app.post("/api/query")
async def query_endpoint(request: QueryRequest):
    """
    HTTP endpoint for queries (non-streaming, returns final result only).
    Use WebSocket for real-time updates.
    """
    try:
        result = await agent.run_streaming(
            request.query,
            session_id=request.session_id
        )
        
        # Extract solution
        import re
        solution_match = re.search(r"<solution>(.*?)</solution>", result.get("solution", ""), re.S)
        solution_text = solution_match.group(1).strip() if solution_match else result.get("solution", "")
        
        return {
            "success": True,
            "session_id": result.get("session_id"),
            "selected_data": result.get("selected_data", {}).get("name", []),
            "selected_tools": result.get("selected_tools", []),
            "solution": solution_text,
            "visualization": result.get("visualization_base64"),
            "visualization_path": result.get("visualization_path")
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/tools")
async def get_available_tools():
    """Get list of available tools."""
    tools_info = {
        name: {
            "description": meta.get("description", ""),
            "args": meta.get("args", {}),
            "returns": meta.get("returns", "")
        }
        for name, meta in TOOLS.items()
    }
    return {"tools": tools_info}


@app.get("/api/data")
async def get_available_data():
    """Get list of available datasets."""
    return {"datasets": DATA_LIST}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=1030, log_level="info")